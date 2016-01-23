#include <iostream>
#include <string>

#include "heat.hpp"
#include "sparse.hpp"

int main(int argc, char *argv[])
{
  /* Get command line arguments */
  if (argc != 3)
  {
    std::cout << "Usage:" << std::endl;
    std::cout << "  " << argv[0] << " <input file> <soln prefix>" << std::endl;
    return 0;
  }
  std::string inputfile   = argv[1];
  std::string soln_prefix   = argv[2];

  // SparseMatrix A;
  // A.AddEntry(0,0,1);
  // A.AddEntry(0,1,2);
  // A.AddEntry(1,1,5);
  // A.AddEntry(2,1,2);
  // A.AddEntry(2,2,1);
  //
  // std::vector<double> a;
  // a.push_back(1);
  // a.push_back(2);
  // a.push_back(3);
  //
  // A.ConvertToCSR();
  // std::vector<double> x = A.MulVec(a);
  //
  // for (int i = 0; i < x.size(); i++) {
  //   std::cout << std::to_string(x[i]) << " ";
  // }
  //   std::cout << std::endl;




  /* Setup 2D heat equation system */
  HeatEquation2D sys;
  int status = sys.Setup(inputfile);
  if (status)
  {
    std::cerr << "ERROR: System setup was unsuccessful!" << std::endl;
    return 1;
  }

  /* Solve system using CG */
  status = sys.Solve(soln_prefix);
  if (status)
  {
    std::cerr << "ERROR: System solve was unsuccessful!" << std::endl;
    return 1;
  }

  return 0;
}
