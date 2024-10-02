/* Kurkov Ivan, 22.B06-MM, 25.09.2024 */
#ifndef CHOL_HPP
#define CHOL_HPP

#include <cmath>

#include "matrix.hpp"
#include "vector.hpp"

/* Calculate Cholesky decomposition for give matrix A
 * @param A - source matrix
 * @param L: LL^T = A
 * @return true: A is positive determined
 * @return false: A isn't positive determined */
template <typename T>
bool CholDecomp( const Matrix<T> &A, Matrix<T> &L )
{
  if (!A.IsSymm())
    return false;

  L = Matrix<T>(A.rows(), A.cols());
  for (size_t i = 0; i < L.rows(); i++)
  {
    T tmp = A[i][i];

    for (size_t j = 0; j < i; j++)
      tmp -= L[i][j] * L[i][j];
    if (tmp < T(0))
      return false;
    L[i][i] = sqrt(tmp);

    for (size_t j = i + 1; j < L.rows(); j++)
    {
      tmp = A[j][i];
      for (size_t k = 0; k < i; k++)
        tmp -= L[j][k] * L[i][k];
      L[j][i] = tmp / L[i][i];
    }
  }
  return true;
}

/* Solve LL^Tx = b
 * @param L - lower triangle matrix from Cholesky decompoition
 * @param b - vector of right sides of system
 * @param opers - number of performed operations
 * @return x - solution of the system */
template <typename T>
Vector<T> CholSolve( const Matrix<T> &L, const Vector<T> &b, size_t &opers )
{
  size_t n = L.rows();
  Vector<T> x(n);

  /* Solve Ly = b */
  for (size_t i = 0; i < n; i++)
  {
    x[i] = b[i];

    for (size_t j = 0; j < i; j++)
      x[i] -= L[i][j] * x[j], opers++;
    x[i] /= L[i][i];
  }
  /* Solve L^Tx = y */
  for (size_t i = n - 1; i + 1 > 0; i--)
  {
    for (size_t j = n - 1; j > i; j--)
      x[i] -= L[j][i] * x[j], opers++;
    x[i] /= L[i][i];
  }
  return x;
}

#endif // !CHOL_HPP
