/* Kurkov Ivan, 22.B06-MM, 16.09.2024 */
#ifndef LUD_HPP

#include <stdexcept>

#include "matrix.hpp"
#include "vector.hpp"

/* Find lower and upper triangle matrices L and U such that A = LU */
/* @return Determinant of A*/
template <typename T>
T LUDecomp( const Matrix<T> &A, Matrix<T> &L, Matrix<T> &U )
{
  size_t n = A.rows();
  T det = T(1);

  L = AlmUnitMatrix<T>(n, n);
  U = Matrix<T>(n, n);

  for (size_t i = 0; i < n; i++)
    for (size_t j = 0; j < n; j++)
    {
      T sum = T(0);

      if (i <= j)
      {
        for (size_t k = 0; k <= i; k++)
          sum += L[i][k] * U[k][j];
        U[i][j] = A[i][j] - sum;
      }
      else
      {
        for (size_t k = 0; k <= j; k++)
          sum += L[i][k] * U[k][j];
        L[i][j] = (A[i][j] - sum) / U[j][j];
      }
    }
  for (size_t i = 0; i < n; i++)
    det *= U[i][i];
  return det;
}

/* Solve LUx = b where L and U - lower and upper trangle matrixes respectively */
/* @return Vector x - solution of the system */
template <typename T>
Vector<T> LUSolve( const Matrix<T> &L, const Matrix<T> &U, const Vector<T> &b, size_t &opers )
{
  size_t n = L.rows();
  Vector<T> x(n);

  opers = 0;
  /* Solve Ly = b */
  for (size_t i = 0; i < n; i++)
  {
    x[i] = b[i];

    for (size_t j = 0; j < i; j++)
      x[i] -= L[i][j] * x[j], opers++;
  }
  /* Solve Ux = y */
  for (size_t i = n - 1; i + 1 > 0; i--)
  {
    for (size_t j = n - 1; j > i; j--)
      x[i] -= U[i][j] * x[j], opers++;
    x[i] /= U[i][i];
  }
  return x;
}

/* Find triangle matrices L, U and permutation matrix P such that PA = LU */
/* @return Determinant of A*/
template <typename T>
T LUPDecomp( const Matrix<T> &A, Matrix<T> &L, Matrix<T> &U, Matrix<T> &P )
{
  size_t n = A.rows(), *perm = new size_t[n], cnt = 0;
  T det = T(1);

  for (size_t i = 0; i < n; i++)
    perm[i] = i;
  L = AlmUnitMatrix<T>(n, n);
  U = A;
  P = Matrix<T>(n, n);

  /* Compute U = U + L - E */
  for (size_t i = 0; i < n; i++)
  {
    T max = abs(U[i][i]);
    size_t idx = i;

    for (size_t j = i + 1; j < n; j++)
      if (abs(U[j][i]) > max)
        max = abs(U[j][i]), idx = j;
    if (max == 0)
      throw std::invalid_argument("Matrix A isn't invertible!");
    std::swap(perm[i], perm[idx]);
    for (size_t j = 0; j < n; j++)
      std::swap(U[i][j], U[idx][j]);
    if (i != idx)
      cnt++;

    for (size_t j = i + 1; j < n; j++)
    {
      U[j][i] /= U[i][i];
      for (size_t k = i + 1; k < n; k++)
        U[j][k] -= U[j][i] * U[i][k];
    }
  }
  /* Separate U into L and U */
  for (size_t i = 1; i < n; i++)
    for (size_t j = 0; j < i; j++)
    {
      L[i][j] = U[i][j];
      U[i][j] = T(0);
    }
  for (size_t i = 0; i < n; i++)
    det *= U[i][i];
  for (size_t i = 0; i < n; i++)
    P[i][perm[i]] = T(1);
  return det * (cnt & 1 ? -1 : 1);
}

/* Solve LUx = Pb where L and U - lower and upper trangle matrixes respectively */
/* @return Vector x - solution of the system */
template <typename T>
Vector<T> LUPSolve( const Matrix<T> &L, const Matrix<T> &U, const Matrix<T> &P, const Vector<T> &b, size_t &opers )
{
  return LUSolve(L, U, P * b, opers);
}

template <typename T>
T GaussElimination( Matrix<T> A, const Vector<T> &b, Vector<T> &x, size_t &opers )
{
  using row_t = typename Matrix<T>::Row;

  row_t *rows = new row_t[A.rows()];
  T det = 1;
  size_t n = std::min(A.rows(), A.cols());

  for (size_t i = 0; i < A.rows(); i++)
    rows[i] = A[i];
  x = b;

  opers = 0;
  /* Forward move */
  for (size_t i = 0; i < n; i++)
  {
    T mx = abs(rows[i][i]), tmp;
    size_t idx = i;

    /* Select pivot element from the column */
    for (size_t j = i + 1; j < A.rows(); j++)
      if (mx < abs(rows[j][i]))
      {
        mx = abs(rows[j][i]);
        idx = j;
      }
    if (mx == 0)
      throw std::invalid_argument("System hasn't got solutions");
    /* Swap rows */
    swap<T>(rows[i], rows[idx]);
    std::swap(x[i], x[idx]);

    tmp = rows[i][i];
    det *= tmp * (i != idx ? -1 : 1);

    for (size_t j = i + 1; j < A.cols(); j++)
      rows[i][j] /= tmp;
    x[i] /= tmp;
    for (size_t j = i + 1; j < A.rows(); j++)
    {
      tmp = rows[j][i];
      for (size_t k = i + 1; k < A.cols(); k++)
        rows[j][k] -= rows[i][k] * tmp, opers++;
      x[j] -= x[i] * tmp;
    }
  }
  /* Backward move */
  for (size_t i = n - 1; i > 0; i--)
    for (size_t j = n - 1; j >= i; j--)
      x[i - 1] -= rows[i - 1][j] * x[j], opers++;
  return det;
}

#endif // !LUD_HPP