pyvoro
======

3D Voronoi tessellations: a python entry point for the [voro++ library](http://math.lbl.gov/voro++/)
by Chris H. Rycroft (UC Berkeley / Lawrence Berkeley Laboratory). This package is based on [pyvoro](https://github.com/joe-jordan/pyvoro) python package
by Joe Jordan (Imperial College London) but it has been completely revamped.


Installation
------------

Installation from source is the same as for any other python module. Issuing 
  
    python setup.py install
    
will install pyvoro system-wide, while 

    python setup.py install --user

will install it only for the current user. Any 
[other](https://pythonhosted.org/an_example_pypi_project/setuptools.html#using-setup-py)  `setup.py` keywords 
can also be used, including 
 
    python setup.py develop
    
to install the package in 'development' mode. Alternatively, if you want all the dependencies pulled in automatically,  
you can still use `pip`:

    pip install -e .

`-e` option makes pip install package from source in development mode. 

You can then use the code with:

    from pyvoro import Voronoi
    Voronoi( ... )

Example:
--------

```python
from pyvoro import Voronoi
voro = Voronoi(
  [[1.0, 2.0, 3.0], [4.0, 5.5, 6.0]],
  [[0.0, 10.0], [0.0, 10.0], [0.0, 10.0]],
  radii=[1.3, 1.4]
)
```

creating a Voronoi object with cells calculated, in particular their volume :

```python
voro.volumes
array([207.10205186, 792.89794814])
```

Cells connectivity with others can be retrieve with the `connectivity` attribute:

```python
cindex, cvalues = voro.connectivity
for start, end in zip(cindex[:-1], cindex[1:]):
      print(cvalues[start:end])
[-5 -2 -3  1 -1 -6]
[ 0 -3 -6 -5 -1 -4 -2]
```

Note that **negative ids** correspond to the boundaries (of which
there are three at the corner of a box, specifically ids `1`, `3` and `5`, (the
`x_i = 0` boundaries, represented with negative ids hence `-1`, `-3` and `-5` --
this is voro++'s conventional way of referring to boundary interfaces.)
