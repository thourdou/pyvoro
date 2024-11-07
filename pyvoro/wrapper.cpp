#include "wrapper.h"
#include "../src/voro++.hh"
#include <stdio.h>
using namespace voro;
using namespace std;

void* container_poly_create(double ax_, double bx_, double ay_, double by_,
  double az_, double bz_, int nx_, int ny_, int nz_, int px_, int py_, int pz_) {
  
  return (void*) new container_poly(ax_, bx_, ay_, by_, az_, bz_, nx_, ny_, nz_, (bool)px_,
      (bool)py_, (bool)pz_, 8);
}


void put_walls(void* con, int nwalls, int* wids, double* wv){
  container_poly* c = (container_poly*)con;
  int v = 0;

  /* Hardcode wall ids from pyx:
  Plane=0, Sphere=1, Cylinder=2, Cone=3
  */
  for (int i = 0; i < nwalls; i++) {
    if (wids[i] == 0){
      c->add_wall(new wall_plane(wv[v],wv[v+1],wv[v+2],wv[v+3],-7-i));
      v=v+4;
    }
    else if (wids[i] == 1){
      c->add_wall(new wall_sphere(wv[v],wv[v+1],wv[v+2],wv[v+3],-7-i));
      v=v+4;
    }
    else if (wids[i] == 2){
      c->add_wall(new wall_cylinder(wv[v],wv[v+1],wv[v+2],wv[v+3],wv[v+4],wv[v+5],wv[v+6],-7-i));
      v=v+7;
    }
    else if (wids[i] == 3){
      c->add_wall(new wall_cone(wv[v],wv[v+1],wv[v+2],wv[v+3],wv[v+4],wv[v+5],wv[v+6],-7-i));
      v=v+7;
    }
    else{
      printf("Unknown wall type\n");
      break;
    }
  }
}


void compute_voronoi(void* container_poly_, int n_, double* points, double* radii,
  double* volumes, double* centers, int* cindex, std::vector<int>* cvalues) {

  container_poly* con = (container_poly*)container_poly_;
  particle_order vo;
  c_loop_order* cla = new c_loop_order(*(con),vo);

  voronoicell_neighbor cell;
  std::vector<int> vi;
  std::vector<double> positions;

  double x, y, z, r;
  int i, v;

  // Fill container with particles 
  for (i = 0; i < n_; i++) {
    con->put(vo, i, points[3*i], points[3*i+1], points[3*i+2], radii[i]);
  }

  // Initialize outputs
  cvalues->clear();
  cindex[0] = 0;

  if(cla->start()) do if (con->compute_cell(cell, *(cla))) {
    // Get information for the current particle
    cla->pos(i, x, y, z, r);

    volumes[i] = cell.volume();

    cell.centroid(centers[3*i], centers[3*i+1], centers[3*i+2]);
    centers[3*i]+=x;
    centers[3*i+1]+=y;
    centers[3*i+2]+=z;

    cell.neighbors(vi);
    cindex[i+1] = vi.size();
    for (v=0; v<vi.size(); v++) cvalues->push_back(vi[v]);
    
  } while (cla->inc());

  // Cumulative sum of connectivity index
  for (i = 0; i < n_; i++) {
    cindex[i+1] += cindex[i];
  }

  delete cla;
}


void find_voronoi_cell(void* container, int npoints, double* points, int* cells_ids){
  container_poly* con = (container_poly*)container;
  int pid;
  double x, y, z;

  for (int i=0; i<npoints; i++){
    if ( con->find_voronoi_cell(points[3*i], points[3*i+1], points[3*i+2], x, y, z, pid) ){
      cells_ids[i] = pid;
    }
    else cells_ids[i] = -1;
  }
}


void dispose_all(void* container_poly_) {
  delete (container_poly*)container_poly_;
}

