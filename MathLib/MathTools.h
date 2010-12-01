/*
 * \file MathTools.h
 *
 *  Created on: Jan 13, 2010
 *      Author: fischeth
 */

#ifndef MATHTOOLS_H_
#define MATHTOOLS_H_

#include <vector>
#include <cmath>
#include <limits>
#include "Point.h"

namespace MATHLIB {

/**
 * standard inner product in R^3
 * \param v0 array of type T representing the vector
 * \param v1 array of type T representing the vector
 * \param n the size of the array
 * */
template<class T>
double scpr(const T* v0, const T* v1, size_t n) {
	long double res(0.0);
	for (size_t k(0); k < n; k++)
		res += v0[k] * v1[k];
	return (double) res;
}

/**
 * computes the cross (or vector) product of the 3d vectors u and v
 * the result is given in the vector r
 */
void crossProd (const double u[3], const double v[3], double r[3]);

/**
 * calcProjPntToLineAndDists computes the orthogonal projection
 * of a point p to the line described by the points a and b,
 * $g(\lambda) = a + \lambda (b - a)$,
 * the distance between p and the projected point
 * and the distances between the projected point and the end
 * points a, b of the line
 * \param p the (mesh) point
 * \param a first point of line
 * \param b second point of line
 * \param lambda the projected point described by the line equation above
 * \param d0 distance to the line point a
 * \returns the distance between p and the orthogonal projection of p
 */
double calcProjPntToLineAndDists(const double p[3], const double a[3],
		const double b[3], double &lambda, double &d0);

/** squared euklid norm of the vector p0 */
double sqrNrm2(const GEOLIB::Point* const p0);

/** squared dist between GEOLIB::Points p0 and p1 */
double sqrDist(const GEOLIB::Point* p0, const GEOLIB::Point* p1);

/** squared dist between double arrays p0 and p1 (size of arrays is 3) */
double sqrDist(const double* p0, const double* p1);

/** linear normalisation of val from [min, max] into [0,1] */
float normalize(const float min, const float max, const float val);

} // namespace

#endif /* MATHTOOLS_H_ */
