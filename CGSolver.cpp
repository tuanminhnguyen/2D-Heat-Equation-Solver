#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include "matvecops.hpp"
#include "sparse.hpp"

/* Write solution to a file
    - soln_prefix is the solution prefix of the solution files
    - nx is number of grid points in x direction
    - sol is a vector storing the solution at any iteration
    - upbnd is a vector storing the hot boundary temperatures
    - lowbnd is a vector storing the cold boundary temperatures
==============================================================================*/
void writeSoln(std::string         soln_prefix,
               uint                 niter,
               int              nx,
               std::vector<double> &sol,
               std::vector<double> &upbnd,
               std::vector<double> &lowbnd,
               double hs) {

   std::string niterStr;
   // add 0's to prefix of the solution file, accordingly
   if (niter < 10) {
     niterStr = "00" + std::to_string(niter);
   }
   else if (niter < 100) {
     niterStr = "0" + std::to_string(niter);
   }
   else {
     niterStr = std::to_string(niter);
   }

   std::ofstream f(soln_prefix + niterStr + ".txt");

   // write the hot boundary temperatures
   // --------------------------------------------------------------------------
   for (int i = 0; i < upbnd.size(); i++) {
     f << upbnd[i] << std::setfill(' ') << std::setw(8);
   }
   f << std::endl;

  //  write the interior nodes' temperatures
  // ---------------------------------------------------------------------------
   for (int i = 0; i < sol.size(); i++) {
     if ( (i+1) % nx == 0 ) {
       f << sol[i] << std::endl;
     }
     else {
       f << std::to_string(sol[i]) << std::setfill(' ') << std::setw(12);
     }
   }

   // write the cold boundary temperatures
   // --------------------------------------------------------------------------
   for (int i = 0; i < lowbnd.size(); i++) {
     f << lowbnd[i] << std::setfill(' ') << std::setw(12);
   }

   f.close();
}


/* Compute an approximate solution to the equation Ax = b using
Conjugate Gradient method.
Inputs:
   - A is a SparseMatrix object
   - b is the right hand side vector, in full form, not in sparse format
   - x is the initial guess for the solution
   - tol is the tolerance level for acceptable solution
   - nx and ny are numbers of grid points in x and y directions
   - hs is the discretization size
==============================================================================*/
int CGSolver(SparseMatrix       &A,
            std::vector<double> &b,
            std::vector<double> &x,
            double              tol,
            std::string         soln_prefix,
            int                 nx,
            int                 ny,
            double              hs) {

// initialize variables: un (storing x_n), pn (p_n and p_{n+1} are not needed
// in the same iteration, so no need to use a pn1 variable),
// rn and rn1 (r_n and r_{n+1}, both are needed in the same iteration)
std::vector<double> un, pn, rn, rn1, upbnd, lowbnd,
                    temp1, temp2, temp3, temp4, temp5;
// These temp's variables are needed to 'anchor' the memory addresses;
// Otherwise compiler complains about l-value required.

double L2normr0 = 0.0, L2normr = 0.0, alpha = 0.0, beta = 0.0;
uint niter = 0, nitermax = (uint) b.size();

for (int i = 0; i < nx; i++) {
    upbnd.push_back(b[i]/hs);
    lowbnd.push_back(b[nx*ny - nx + i]/hs);
}

un = x; //un now stores the initial guess
temp1 = A.MulVec(x);
rn = vectorSubtract(b, temp1);
L2normr0 = twoNorm(rn);
pn = rn;

// conjugate gradient algorithm
// -----------------------------------------------------------------------------
while (niter < nitermax) {
  niter++;
  temp2 = A.MulVec(pn);
  alpha = dotProduct(rn,rn) / dotProduct(pn, temp2);

  // after execution of this line, x stores the most recently updated solution
  temp3 = scalarMult(alpha, pn);
  x = vectorAdd(un, temp3);
  temp2 = A.MulVec(pn);
  temp4 = scalarMult(alpha, temp2);
  rn1 = vectorSubtract(rn, temp4);
  L2normr = twoNorm(rn1);

  if ( niter % 10 == 0 ) {
    writeSoln(soln_prefix, niter, nx, x, upbnd, lowbnd, hs);
  }

  // if the norm ratio satisfies the tolerance, terminate conjugate gradient
  // and write the final solution to file
  if ((L2normr / L2normr0) < tol) {
    std::cout << "SUCCESS: CG solver converged in " << niter << " iterations." << std::endl;
    writeSoln(soln_prefix, niter, nx, x, upbnd, lowbnd, hs);
    return 0;
  }

  beta = dotProduct(rn1,rn1) / dotProduct(rn,rn);
  temp5 = scalarMult(beta,pn);
  pn = vectorAdd(rn1, temp5);

  un = x;
  rn = rn1;
  }

// If niter == nitermax - 1, CG did not converge, so return -1
return -1;
}
