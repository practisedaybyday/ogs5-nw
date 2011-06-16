/*
 * TriangularSolve.h
 *
 *  Created on: May 6, 2010
 *      Author: TF
 */

#ifndef TRIANGULARSOLVE_H_
#define TRIANGULARSOLVE_H_

namespace MathLib {

/**
 * solves the \f$n \times n\f$ triangular linear system \f$L \cdot y = b\f$,
 * assumes \f$L_{ii} = 1.0\f$, \f$i=1,...,n\f$, \f$b\f$ is destroyed
 * @param L the lower triangular matrix
 * @param b at beginning the right hand side vector, at the end the solution vector
 */
void forwardSolve (const Matrix <double> &L, double* b);

/**
 * solves the \f$n \times n\f$ triangular linear system \f$U \cdot x=y\f$,
 * \f$U\f$, where \f$U\f$ is a upper triangular matrix.
 * @param U upper triangular matrix
 * @param y at beginning the right hand side, at the end the solution
 */
void backwardSolve (const Matrix <double> &U, double* y);

} // end namespace MathLib

#endif /* TRIANGULARSOLVE_H_ */
