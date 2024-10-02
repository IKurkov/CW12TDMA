/* Kurkov Ivan, 22.B06-MM, 02.10.2024 */
#include <conio.h>
#include <cstdlib>
#include <iostream>

#include "chol.hpp"
#include "fort.hpp"
#include "lud.hpp"
#include "matrix.hpp"
#include "vector.hpp"

#include "tdma.h"

Matrix<double> GenThreeDiag( size_t n )
{
  Matrix<double> A(n, n);

  A[0][0] = A[0][1] = (10.0 * rand() / RAND_MAX + 10.0) * (rand() & 1 ? 1.0 : -1.0);
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
  Matrix<double> A = GenThreeDiag(3);
  Vector<double> b(3), x;

  b[1] = 2;
  std::cout << A << '\n';
  if (ThreeDiagMatrixAlg(A, b, x))
  {
    std::cout << "Ok\n" << NormInf(A * x - b);
  }
  else
    std::cout << "Bad arguments\n";

  return 0;
}