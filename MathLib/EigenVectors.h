/**
 * \file EigenVector.h
 * 11/01/2012 LB Initial implementation
 * Shamelessly stolen from VTK: Common/vtkMath.h / .cxx
 */

#ifndef EIGENVECTOR_H
#define EIGENVECTOR_H

namespace MathLib
{

/**
 * Jacobi iteration for the solution of eigenvectors/eigenvalues of a nxn
 * real symmetric matrix. Square nxn matrix a; size of matrix in n;
 * output eigenvalues in w; and output eigenvectors in v. Resulting
 * eigenvalues/vectors are sorted in decreasing order; eigenvectors are
**/
template<class T> int JacobiN(T **a, int n, T *w, T **v);

int JacobiN(float **a, int n, float *w, float **v);
int JacobiN(double **a, int n, double *w, double **v);

/**
 * Jacobi iteration for the solution of eigenvectors/eigenvalues of a 3x3
 * real symmetric matrix. Square 3x3 matrix a; output eigenvalues in w;
 * and output eigenvectors in v. Resulting eigenvalues/vectors are sorted
 * in decreasing order; eigenvectors are normalized.
**/
int Jacobi(float **a, float *w, float **v);
int Jacobi(double **a, double *w, double **v);

} // namespace MathLib

#endif // EIGENVECTOR_H
