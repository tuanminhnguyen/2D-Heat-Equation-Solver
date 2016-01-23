#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include "matvecops.hpp"

/* Multiply a vector by a scalar */
std::vector<double> scalarMult(double &val, std::vector<double> &vec) {
  std::vector<double> ans;
  for (uint i = 0; i < vec.size(); ++i) {
    ans.push_back(vec[i] * val);
  }
  return ans;
  }

/* Add two vectors of the same size */
std::vector<double> vectorAdd(std::vector<double> &x, std::vector<double> &y) {
  std::vector<double> ans;
  for (uint i = 0; i < x.size(); ++i) {
    ans.push_back(x[i] + y[i]);
  }
  return ans;
  }

/* Subtract two vectors of the same size */
std::vector<double> vectorSubtract(std::vector<double> &x, std::vector<double> &y) {
  std::vector<double> ans;
  for (uint i = 0; i < x.size(); ++i) {
    ans.push_back(x[i] - y[i]);
  }
  return ans;
  }

/* Compute the inner product of two vectors of the same size */
double dotProduct(std::vector<double> &x, std::vector<double> &y) {
  double ans = 0.0;
  for (uint i = 0; i < x.size(); ++i) {
   ans += (x[i] * y[i]);
  }
  return ans;
  }

/* Compute the L2-norm of a vector */
double twoNorm(std::vector<double> &x) {
    return sqrt(dotProduct(x,x));
}
