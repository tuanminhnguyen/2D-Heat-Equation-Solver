#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include "heat.hpp"
#include "CGSolver.hpp"
#include "sparse.hpp"

/* Method to setup Ax=b system
==============================================================================*/
int HeatEquation2D::Setup(std::string inputfile) {
  std::ifstream input;
  input.open(inputfile.c_str());
  double L, W, h, Tc, Th;

  if ( !(input >> L >> W >> h >> Tc >> Th) ) {
    std::cout << "Cannot open " << inputfile << std::endl;
    return -1;
  }

  nx = (int) ceil(L / h);     // # distinct gridpts in x
  ny = (int) ceil(W / h) - 1; // # distinct gridpts in y, excluding boundaries
  int t = nx * ny;
  hs = 1/(h*h);
  double a = -1.0*hs;
  double c = 4.0*hs;


  // Set up vector x
  // ---------------------------------------------------------------------------
  x.resize((size_t) t, 1.0); // nx * ny is number of unknowns

  // Set up vector b
  // ---------------------------------------------------------------------------
  b.resize((size_t) t, 0.0);
  // Upper boundary, constant temperature
  for (uint i = 0; i < (uint) nx; i++) { b[i] = (-1.0) * a * Th; }

  // Lower boundary, temperature follows given function
  double x = 0;
  for (uint i = (uint) (t - nx); i < (uint) t; i++) {
    b[i] = a*Tc*(exp(-10.0 * pow(x-L/2, 2)) -2);
    x += h;
  }

  // Set up matrix A
  // First block, from row 0 to row nx-1
  // ---------------------------------------------------------------------------
  A.AddEntry(0, 0,    c);
  A.AddEntry(0, 1,    a);
  A.AddEntry(0, nx-1, a);
  A.AddEntry(0, nx,   a);

  for (int i = 1; i < nx-1; i++) {
    A.AddEntry(i, i-1,  a);
    A.AddEntry(i, i,    c);
    A.AddEntry(i, i+1,  a);
    A.AddEntry(i, i+nx, a);
  }

  A.AddEntry(nx-1, 0,      a);
  A.AddEntry(nx-1, nx-2,   a);
  A.AddEntry(nx-1, nx-1,   c);
  A.AddEntry(nx-1, 2*nx-1, a);

  // Second block to last to second block
  // ---------------------------------------------------------------------------
  for (int i = nx; i < t-nx; i++) {
    // A.AddEntry(i, i-nx, a);
    if (i%nx == nx-1) {
      A.AddEntry(i, i-nx, a);
      A.AddEntry(i, i-nx+1, a);
      A.AddEntry(i, i-1,    a);
      A.AddEntry(i, i, c);
      A.AddEntry(i, i+nx,a );
    }
    else if (i%nx == 0) {
      A.AddEntry(i, i-nx, a);
      A.AddEntry(i, i, c);
      A.AddEntry(i, i+1,    a);
      A.AddEntry(i, i+nx-1, a);
      A.AddEntry(i, i+nx,a );}
    else {
      A.AddEntry(i, i-nx, a);
      A.AddEntry(i, i-1, a);
      A.AddEntry(i, i, c);
      A.AddEntry(i, i+1,    a);
      A.AddEntry(i, i+nx,a );
    }
  }

  // Last block
  // ---------------------------------------------------------------------------
  A.AddEntry(t-nx, t-2*nx, a);
  A.AddEntry(t-nx, t-nx,   c);
  A.AddEntry(t-nx, t-nx+1, a);
  A.AddEntry(t-nx, t-1,      a);

  for (int i = t-nx+1; i < t-1; i++) {
    A.AddEntry(i, i-nx, a);
    A.AddEntry(i, i-1,  a);
    A.AddEntry(i, i,    c);
    A.AddEntry(i, i+1,  a);
  }

  A.AddEntry(t-1, t-nx-1, a);
  A.AddEntry(t-1, t-nx,   a);
  A.AddEntry(t-1, t-2,    a);
  A.AddEntry(t-1, t-1,    c);
  A.ConvertToCSR();
  // A.tprint();
  return 0;
}

/* Method to solve system using CGsolver
==============================================================================*/
int HeatEquation2D::Solve(std::string soln_prefix) {
   CGSolver(A, b, x, 1.0e-5, soln_prefix, nx, ny, hs);
   return 0;
}
