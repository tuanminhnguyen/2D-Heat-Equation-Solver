#ifndef SPARSE_HPP
#define SPARSE_HPP

#include <vector>
#include <iostream>
#include "COO2CSR.hpp"

class SparseMatrix
{
  private:
    std::vector<int> i_idx;
    std::vector<int> j_idx;
    std::vector<double> a;    


  public:
    /* Method to add entry to matrix in COO format */
    void AddEntry(int i, int j, double val);

    /* Method to convert COO matrix to CSR format using provided function */
    void ConvertToCSR();

    /* Method to perform sparse matrix vector multiplication using CSR formatted matrix */
    std::vector<double> MulVec(std::vector<double> &vec);

};

#endif /* SPARSE_HPP */
