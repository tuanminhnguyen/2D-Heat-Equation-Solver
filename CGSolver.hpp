#ifndef CGSOLVER_HPP
#define CGSOLVER_HPP

#include <vector>
#include "sparse.hpp"

/* Function that implements the CG algorithm for a linear system

   Ax = b

   where A is in CSR format.  The starting guess for the solution
   is provided in x, and the solver runs a maximum number of iterations
   equal to the size of the linear system.  Function returns the
   number of iterations to converge the solution to the specified
   tolerance, or -1 if the solver did not converge. */

void writeSoln(std::string         soln_prefix,
             uint                 niter,
             int              nx,
             std::vector<double> &sol,
             std::vector<double> &upbnd,
             std::vector<double> &lowbnd,
             double hs);

int CGSolver(SparseMatrix        &A,
             std::vector<double> &b,
             std::vector<double> &x,
             double              tol,
             std::string soln_prefix,
             int nx, int ny, double hs);

#endif /* CGSOLVER_HPP */
