/* Kurkov Ivan, 22.B06-MM, 02.10.2024 */
#include <conio.h>
#include <cstdlib>
#include <iostream>

#include "fort.hpp"
#include "lud.hpp"
#include "matrix.hpp"
#include "vector.hpp"

#include "tdma.h"

Matrix<double> GenThreeDiag( size_t n )
{
  Matrix<double> A(n, n);

  A[0][0] = (10.0 * rand() / RAND_MAX + 10.0) * (rand() & 1 ? 1.0 : -1.0);
  A[0][1] = 10.0 * rand() / RAND_MAX - 5;
  for (size_t i = 1; i + 1 < n; i++)
  {
    A[i][i - 1] = 10 * rand() / RAND_MAX - 5;
    A[i][i] = (10.0 * rand() / RAND_MAX + 10.0) * (rand() & 1 ? 1.0 : -1.0);
    A[i][i + 1] = 10 * rand() / RAND_MAX - 5;
  }
  A[n - 1][n - 2] = 10.0 * rand() / RAND_MAX - 5;
  A[n - 1][n - 1] = (10.0 * rand() / RAND_MAX + 10.0) * (rand() & 1 ? 1.0 : -1.0);
  return A;
}

int main( void )
{
  bool run = true;

  while (run)
  {
    std::cout << "=====Threediagonal matrix algorithm menu=====\n"
      << "0 - exit\n"
      << "1 - threediagonal matrix algorithm\n";
    switch (_getch())
    {
    case '0':
      run = false;
      break;
    case '1':
    {
      size_t n, ops;
      Matrix<double> A, L, U, P;
      Vector<double> b, x, x_ref;
      fort::char_table compare;

      std::cout << "Input size of matrix: ";
      std::cin >> n;
      A = GenThreeDiag(n);
      std::cout << "A = " << A << "\n";

      x_ref = Vector<double>(n);
      for (size_t i = 0; i < n; i++)
        x_ref[i] = 1;
      std::cout << "Referens solution: " << x_ref << "\n";
      b = A * x_ref;

      compare << fort::header << "Method" << "x" << "|x - x_ref|_inf" << "Number of operations" << fort::endr;

      GaussElimination(A, b, x, ops);
      compare << "Gauss" << x << NormInf(x - x_ref) << ops << fort::endr;

      LUPDecomp(A, L, U, P);
      LUPSolve(L, U, P, b, ops);
      compare << "LUP" << x << NormInf(x - x_ref) << ops << fort::endr;

      ThreeDiagMatrixAlg(A, b, x, ops);
      compare << "Threediagonal" << x << NormInf(x - x_ref) << ops << fort::endr;
      std::cout << compare.to_string();
      break;
    }
    default:
      std::cout << "[Error]: Incorrect choice!\n";
    }
  }
  return 0;
}