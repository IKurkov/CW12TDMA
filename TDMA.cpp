/* Kurkov Ivan, 22.B06-MM, 02.10.2024 */
#include "matrix.hpp"
#include "vector.hpp"

#include "tdma.h"

bool IsDiagDomin( const Matrix<double> &A )
{
  size_t n = A.rows();

  if (fabs(A[0][0]) < fabs(A[0][1]))
    return false;
  for (size_t i = 1; i + 1 < n; i++)
    if (fabs(A[i][i]) < fabs(A[i][i - 1]) + fabs(A[i][i + 1]) || fabs(A[i][i]) <= fabs(A[i][i - 1]))
      return false;
  if (fabs(A[n - 1][n - 1]) <= fabs(A[n - 1][n - 2]))
    return false;

  return true;
}

bool ThreeDiagMatrixAlg( const Matrix<double> &A, const Vector<double> b, Vector<double> &x, size_t &opers )
{
  if (!IsDiagDomin(A))
    return false;

  size_t n = A.rows();
  double *alpha = new double[n - 1],
         *beta  = new double[n - 1], gamma;

  x = Vector<double>(n);

  opers = 0;
  gamma = A[0][0];
  beta[0] = b[0] / gamma;
  alpha[0] = -A[0][1] / gamma;
  opers += 2;

  /* Forward move */
  for (size_t i = 1; i + 1 < n; i++)
  {
    gamma = A[i][i] + A[i][i - 1] * alpha[i - 1];
    beta[i] = (b[i] - A[i][i - 1] * beta[i - 1]) / gamma;
    alpha[i] = -A[i][i + 1] / gamma;
    opers += 4;
  }
  /* Backward move */
  x[n - 1] = (b[n - 1] - A[n - 1][n - 2] * beta[n - 2])
    / (A[n - 1][n - 1] + A[n - 1][n - 2] * alpha[n - 2]);
  opers += 3;
  for (size_t i = n - 2; i + 1 > 0; i--)
    x[i] = alpha[i] * x[i + 1] + beta[i], opers++;

  delete[] alpha;
  delete[] beta;

  return true;
}