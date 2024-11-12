from libcpp cimport bool
from libcpp.vector cimport vector

import numpy as _np
cimport numpy as _np
_np.import_array()

cdef extern from "wrapper.h":
    void* container_poly_create(double ax_, double bx_, double ay_, double by_,
        double az_, double bz_, int nx_, int ny_, int nz_, bool px_, bool py_, bool pz_)
    void put_walls(void* con, int nwalls, int* wids, double* wvalues)
    void compute_voronoi(void* container_poly_, int n_ , double* points, double* r,
        double* volumes, double* centers, int* cindex, vector[int]* cvalues)

    void find_voronoi_cell(void* container, int npoints, double* points, int* cells_ids)

    void dispose_all(void* container_poly_)



def _get_vector(vec, normalize=False):
    '''
    Internal use. Convert iterable vector to numpy array
    '''
    vec = _np.asarray(vec, dtype=float).squeeze()
    if vec.shape not in ((2,), (3,)):
        raise ValueError('Wrong vector dimension')
    if normalize:
        return vec / _np.linalg.norm(vec)
    return vec


class Plane:
    '''
    Plane wall defined by its normal vector and distance from origin

    Attributes
    ----------
    normal: numpy array
        Normal vector to the plane, pointing outward

    distance: float
        Distance of the plane from the origin

    values: tuple(float)
        Numerical values representing the plane
        (normal[0], normal[1], normal[2], distance)
    '''
    def __init__(self, normal, distance):
        self.normal = _get_vector(normal, normalize=True)
        self.distance = float(distance)
        self.values = (*self.normal, self.distance)


class Sphere:
    '''
    Sphere wall defined by its position vector and radius

    Attributes
    ----------
    position: numpy array
        Position vector of the sphere center

    radius: float
        Radius of the sphere

    values: tuple(float)
        Numerical values representing the sphere
        (position[0], position[1], position[2], radius)
    '''
    def __init__(self, position, radius):
        self.position = _get_vector(position)
        self.radius = float(radius)
        self.values = (*self.position, self.radius)


class Cylinder:
    '''
    Cylinder wall defined by its axis position, axis direction and radius

    Attributes
    ----------
    position: numpy array
        Point on the cylinder central axis

    direction: numpy array
        Cylinder central axis

    radius: float
        Radius of the cylinder circular section

    values: tuple(float)
        Numerical values representing the cylinder
        (pos[0], pos[1], pos[2], dir[0], dir[1], dir[2], radius)
    '''
    def __init__(self, position, direction, radius):
        self.position = _get_vector(position)
        self.direction = _get_vector(direction, normalize=True)
        self.radius = float(radius)
        self.values = (*self.position, *self.direction, self.radius)


class Cone:
    '''
    Cone wall defined by its apex position, axis direction and half angle

    Attributes
    ----------
    position: numpy array
        Apex position vector

    direction: numpy array
        Cone central axis

    angle: float
        Half angle of the cone apex

    values: tuple(float)
        Numerical values representing the cone
        (pos[0], pos[1], pos[2], dir[0], dir[1], dir[2], angle)
    '''
    def __init__(self, position, direction, angle):
        self.position = _get_vector(position)
        self.direction = _get_vector(direction, normalize=True)
        self.angle = float(angle)
        self.values = (*self.position, *self.direction, self.radius)



cdef class Voronoi:
    cdef readonly int npoints, ndim
    cdef readonly _np.ndarray bounding_box, points, radii
    cdef readonly _np.ndarray walls, periodicity
    cdef void* con

    cdef vector[int] cvalues
    cdef readonly _np.ndarray centroids, volumes
    cdef readonly tuple connectivity

    """
    Voronoi diagrams in 3 dimensions.

    Parameters:
    -----------
    points: array-like of floats, shape (npoints, 3)
        Coordinates of points to construct a Voronoi diagram from

    bounding_box: array-like of floats, shape (3, 2)
        Start and end sizes of the boundary box

    walls: array-like, optional
        List of pyvoro geometries (Plane, Sphere, Cylinder, Cone) defining cutting walls

    radii: array-like of floats, shape (npoints, ), optional
        Spheres radii of points to construct a Radical Voronoi diagram

    periodic: bool or array-like of bool, shape (3, )
        Periodicity of the boundary box along x, y and z axes

    Attributes
    ----------
    npoints: int
        Number of input points

    ndim: int = 3
        Number of dimensions

    bounding_box: ndarray of double, shape (3, 2)
        Start and end sizes of the boundary box

    walls: ndarray of walls
        Cutting walls list

    periodicity: ndarray of bool, shape (3, )
        Periodicity of the boundary box along x, y and z axes

    points: ndarray of double, shape (npoints, ndim)
        Coordinates of input points

    radii: ndarray of double, shape (npoints, )
        Radii of input points

    centroids: ndarray of double, shape (npoints, 3)
        Coordinates of voronoi cells centroids

    volumes: ndarray of double, shape (npoints, )
        Voronoi cells volumes

    connectivity: tuple of ndarrays of int
        Voronoi cells connectivity to neighbours ones.

        Since cell geometry is not specific, number of neighbours is not consistent.
        Connectivity is constructed with an indexer for each cell.

        >>> voro = Voronoi( ... )
        >>> index, neighbours = voro.connectivity
        >>> for i in range(voro.npoints):
        >>>     neighbours[index[i]:index[i+1]]

    """
    def __init__(self, points, bounding_box, walls=None, radii=None, periodic=False):
        ''''''
        # Check input points
        cdef _np.ndarray[_np.float64_t, ndim=2] pts
        pts = _np.ascontiguousarray(points, dtype=_np.float64)
        if not _np.all(_np.isfinite(pts)):
            raise ValueError('Wrong input points values (nan/infinite)')
        self.points = pts
        del points

        self.npoints = self.points.shape[0]
        self.ndim = self.points.shape[1]

        # Check bounding box
        self.bounding_box = _np.ascontiguousarray(bounding_box, dtype=_np.float64)
        if self.bounding_box.shape[0] != self.ndim and self.bounding_box.shape[1] != 2:
            raise ValueError('Wrong bounding box dimension')
        del bounding_box

        # Check walls
        walltypes = (Plane, Sphere, Cylinder, Cone)
        if walls is None:
            walls = []
        for wall in walls:
            if not isinstance(wall, walltypes):
                raise TypeError('Wrong wall type')
        self.walls = _np.array(walls, dtype=object)
        del walls

        # Check radii array
        # TODO: if not radius specified, evaluate blocks
        # with cpp pre_container_poly->guess_optimal.
        # else directly evaluate radius average ?
        # What if distance between spheres really greater than spheres radius ?
        cdef _np.ndarray[_np.float64_t, ndim=1] rad
        if radii is None:
            raise NotImplementedError()
        else:
            rad = _np.ascontiguousarray(radii, dtype=_np.float64)
            if len(rad) != self.npoints:
                raise ValueError('Wrong number of radii values')
            if not _np.all(_np.isfinite(rad)):
                raise ValueError('Wrong input radii values (nan/infinite)')
            self.radii = rad
            del radii

            # Compute number of blocks (Copy of cmd_line.cc l:377)
            blocks = ((
                self.bounding_box[:,1]-self.bounding_box[:,0]
              ) * 0.3/self.radii.mean() + 1).astype(int)

        # Check periodicity
        if periodic in (False, True):
            self.periodicity = _np.array([periodic]*3, dtype=_np.bool_)
        else:
            # Unpack and cast. Will raise an error if not an iterable
            # User can specify longer arrays ... but not taken into account
            self.periodicity = _np.array(
                [periodic[0],periodic[1],periodic[2]],
                dtype=_np.bool_
            )
        del periodic

        # Build the container object
        self.con = <void*> container_poly_create(
            self.bounding_box[0][0],
            self.bounding_box[0][1],
            self.bounding_box[1][0],
            self.bounding_box[1][1],
            self.bounding_box[2][0],
            self.bounding_box[2][1],
            blocks[0],
            blocks[1],
            blocks[2],
            self.periodicity[0],
            self.periodicity[1],
            self.periodicity[2]
        )

        # Add walls to the container
        cdef _np.ndarray[_np.int32_t, ndim=1] wids
        cdef _np.ndarray[_np.float64_t, ndim=1] wval
        if len(self.walls):
            wids = _np.array([
                walltypes.index(wall.__class__) for wall in self.walls
            ], dtype=_np.int32)
            wval = _np.concatenate([wall.values for wall in self.walls], dtype=_np.float64)
            put_walls(self.con, len(self.walls), <int*>wids.data, <double*>wval.data)


        # Create output fixed-size arrays
        cdef _np.ndarray[_np.float64_t, ndim=1] vols # Cells volumes
        vols = _np.empty(self.npoints, dtype=_np.float64)*_np.nan
        cdef _np.ndarray[_np.float64_t, ndim=2] centers # Cells centroids
        centers = _np.empty( (self.npoints, self.ndim), dtype=_np.float64)*_np.nan

        # Create output indexers
        # Points are not computed yet so outside ones are not known.
        # Connectivities must set to zeros at first so cpp compute function
        # can do a cumulative sum of non empty cells connectivities
        cdef _np.ndarray[_np.int32_t, ndim=1] cindex # Neighbourhood connectivity index
        cindex = _np.zeros(self.npoints+1, dtype=_np.int32)

        # Create output unknown-size vectors -> will be cast to numpy after

        # Compute the tessellation
        compute_voronoi(
            self.con, self.npoints,
            <double*>pts.data, <double*>rad.data,
            <double*>vols.data, <double*>centers.data,
            <int*>cindex.data, &self.cvalues,
        )

        # Set class attributes
        self.volumes = vols
        self.centroids = centers

        cdef int[::1] cvals_view = <int[:self.cvalues.size()]> &self.cvalues[0]
        self.connectivity = (cindex, _np.asarray(cvals_view))


    def __dealloc__(self):
        ''''''
        dispose_all(self.con)


    def find_cells(self, points):
        '''
        Find voronoi cells containing the given points

        Parameters
        ----------
        points: array-like
            Points to locate

        Returns
        -------
        cids: ndarray of int
            Indices of cells containing each point.
            Points outside the triangulation get the value -1.
        '''
        cdef _np.ndarray[_np.float64_t, ndim=2] pts
        pts = _np.ascontiguousarray(points, dtype=_np.float64)
        if pts.shape[1] != self.ndim:
            raise ValueError('Wrong input points array dimension')

        cdef _np.ndarray[_np.int32_t, ndim=1] cids
        cids = _np.empty(len(pts), dtype=_np.int32)
        find_voronoi_cell(self.con, len(pts), <double*>pts.data, <int*>cids.data)

        return cids
