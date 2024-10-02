/* Kurkov Ivan, 22.B06-MM, 02.10.2024 */
#ifndef TDMA_H
#define TDMA_H

#include "matrix.hpp"
#include "vector.hpp"

/* Solve system Ax = b using threediagonal matrix algorithm
 * @param A threediagonal matrix
 * @param b vector of right sides of system
 * @param x solution of system
 * @return true - A is diagonal dominant, false - A isn't diagonal dominant */
bool ThreeDiagMatrixAlg( const Matrix<double> &A, const Vector<double> b, Vector<double> &x );

/* Check if given matrix is diagonal dominant
* @param A threediagonal matrix
* @return true - A is diagonal dominant, false - A isn't diagonal dominant */
bool IsDiagDomin( const Matrix<double> &A );

#endif // !TDMA_H

