#ifndef MATVECOPS_HPP
#define MATVECOPS_HPP

#include <vector>

std::vector<double> scalarMult(double &val, std::vector<double> &vec);

std::vector<double> vectorAdd(std::vector<double> &x, std::vector<double> &y);

std::vector<double> vectorSubtract(std::vector<double> &x, std::vector<double> &y);

double dotProduct(std::vector<double> &x, std::vector<double> &y);

double twoNorm(std::vector<double> &x);

#endif /* MATVECOPS_HPP */
