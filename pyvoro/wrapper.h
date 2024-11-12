#include <vector>

void* container_poly_create(double ax_, double bx_, double ay_, double by_,
  double az_, double bz_, int nx_, int ny_, int nz_, bool px_, bool py_, bool pz_);

void put_walls(void* con, int nwalls, int* wids, double* wv);

void compute_voronoi(void* container_poly_, int n_, double* points, double* radii,
  double* volumes, double* centers, int* cindex, std::vector<int>* cvalues);

void find_voronoi_cell(void* container, int npoints, double* points, int* cells_ids);

void dispose_all(void* container_poly_);
