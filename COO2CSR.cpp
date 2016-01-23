#include <algorithm>
#include <vector>

#include "COO2CSR.hpp"

/* Bubble sorting function for sparse matrix format conversion,
   used to sort the entries in one row of the matrix. */

void SortRow(std::vector<int>    &row_idx,
             std::vector<int>    &col_idx,
             std::vector<double> &a,
             int                 start,
             int                 end,
             int                 nz,
             int                 current_row)
{
  for (int i = end - 1; i > start; i--)
  {
    for(int j = start; j < i; j++)
    {
      if (col_idx[j] > col_idx[j+1])
      {
        /* Swap the value and the column index */
        double dt = a[j];
        a[j] = a[j+1];
        a[j+1] = dt;

        int it = col_idx[j];
        col_idx[j] = col_idx[j+1];
        col_idx[j+1] = it;
      }
    }
  }

  /* Accumulate duplicate values and adjust vectors */
  for (int j = start; j<end-1; j++)
  {
    if (col_idx[j] == col_idx[j+1])
    {
      a[j] += a[j+1];

      for (int i = j+1; i<nz-1; i++)
      {
        a[i] = a[i+1];
        col_idx[i] = col_idx[i+1];
      }

      for (unsigned int i = current_row + 1; i<row_idx.size(); i++)
      {
        row_idx[i]--;
      }

      a[nz-1] = 0;
      end--;
      j--;
    }

  }
}

/* In place conversion of square matrix from COO to CSR format */

void COO2CSR(std::vector<double> &val,
             std::vector<int>    &i_idx,
             std::vector<int>    &j_idx)
{
  int nrows = *std::max_element(i_idx.begin(), i_idx.end()) + 1;
  int nnonzeros = (int)val.size();

  int *row_start = new int[nrows+1];
  for (int i = 0; i <= nrows; i++)
  {
    row_start[i] = 0;
  }

  /* Determine row lengths */

  for (int i = 0; i < nnonzeros; i++) row_start[i_idx[i]+1]++;
  for (int i = 0; i < nrows; i++) row_start[i+1] += row_start[i];

  for (int init = 0; init < nnonzeros; )
  {
    double dt = val[init];
    int i = i_idx[init];
    int j = j_idx[init];
    i_idx[init] = -1;

    while (true)
    {
      int i_pos = row_start[i];
      double val_next = val[i_pos];
      int i_next = i_idx[i_pos];
      int j_next = j_idx[i_pos];

      val[i_pos] = dt;
      j_idx[i_pos] = j;
      i_idx[i_pos] = -1;
      row_start[i]++;

      if (i_next < 0) break;

      dt = val_next;
      i = i_next;
      j = j_next;
    }
    init++;

    while ((init < nnonzeros) and (i_idx[init] < 0))
    {
      init++;
    }
  }

  /* Copy row pointer */

  for (int i = 0; i < nrows; i++)
  {
    i_idx[i+1] = row_start[i];
  }
  i_idx[0] = 0;

  /* Sort each row */

  for (int i = 0; i < nrows; i++)
  {
    SortRow(i_idx, j_idx, val, i_idx[i], i_idx[i+1], nnonzeros, i);
  }

  i_idx.resize(nrows+1);

  nnonzeros = i_idx[nrows];
  val.resize(nnonzeros);
  j_idx.resize(nnonzeros);

  delete[] row_start;
}

// #ifdef DEBUG
// #include <iostream>
// int main()
// {
//   std::vector<double> val = {4, 1, 1, 4, 1, 1, 4};
//   std::vector<int> i_idx = {0, 0, 1,  1, 1, 2, 2};
//   std::vector<int> j_idx = {0, 1, 0,  1, 2, 1, 2};
//
//   std::vector<double> val2 = {4, 1, 1, 2, 2, 1, 1, 4};
//   std::vector<int> i_idx2 = {0, 0, 1, 1, 1, 1, 2, 2};
//   std::vector<int> j_idx2 = {0, 1, 0, 1, 1, 2, 1, 2};
//
//
//   COO2CSR(val, i_idx, j_idx);
//   COO2CSR(val2, i_idx2, j_idx2);
//
//
//   /* Matrix without accumulated value */
//   for (auto v : val)
//     std::cout << v << " ";
//   std::cout << std::endl;
//
//   for (auto j : j_idx)
//     std::cout << j << " ";
//   std::cout << std::endl;
//
//   for (auto i : i_idx)
//     std::cout << i << " ";
//   std::cout << std::endl;
//
//   /* Matrix with accumulated value */
//   for (auto v : val2)
//     std::cout << v << " ";
//   std::cout << std::endl;
//
//   for (auto j : j_idx2)
//     std::cout << j << " ";
//   std::cout << std::endl;
//
//   for (auto i : i_idx2)
//     std::cout << i << " ";
//   std::cout << std::endl;
//
//   return 0;
// }
// #endif
