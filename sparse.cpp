#include <vector>
#include "COO2CSR.hpp"
#include "sparse.hpp"

/* Method to add entry to matrix in COO format */
void SparseMatrix::AddEntry(int i, int j, double val){
  i_idx.push_back(i);
  j_idx.push_back(j);
  a.push_back(val);
}

/* Method to convert COO matrix to CSR format using provided function */
void SparseMatrix::ConvertToCSR() {
  COO2CSR(a, i_idx, j_idx);
}


/* Method to perform sparse matrix vector multiplication using CSR formatted
   matrix */
std::vector<double> SparseMatrix::MulVec(std::vector<double> &vec) {
  std::vector<double> ans;

  unsigned int start = 0, end = 0;
  double sum;
  for (unsigned int r = 0; r < i_idx.size() - 1; r++) {
    sum = 0;
    start = (unsigned int) i_idx[r];
    end = (unsigned int) i_idx[r+1];

    for (unsigned int i = start; i < end; i++) {
      if (end > a.size()) {
      end = (unsigned int) a.size();
      }
      sum += a[i] * vec[(unsigned int) j_idx[i]];
    }
    ans.push_back(sum);
  }
  return ans;
}
