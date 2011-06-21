/*
 * \file Triangle.cpp
 *
 *  Created on: Jun 6, 2011
 *      Author: TF
 */

#include "Triangle.h"

// MathLib
#include "LinAlg/GaussAlgorithm.h"
//#include "LinAlg/HouseholderAlgorithm.h"
//#include "LinAlg/TriangularSolve.h"
#include "MathTools.h"
#include "Matrix.h"
#include "Vector3.h"

namespace GEOLIB {

Triangle::Triangle (std::vector<Point *> const &pnt_vec) :
	_pnts(pnt_vec), _initialized (false), _longest_edge (0.0)
{}

Triangle::Triangle (std::vector<Point *> const &pnt_vec, size_t pnt_a, size_t pnt_b, size_t pnt_c) :
	_pnts(pnt_vec), _initialized (true), _longest_edge (0.0)
{
	_pnt_ids[0] = pnt_a;
	_pnt_ids[1] = pnt_b;
	_pnt_ids[2] = pnt_c;
	_longest_edge = MathLib::sqrDist (_pnts[_pnt_ids[0]], _pnts[_pnt_ids[1]]);
	double tmp (MathLib::sqrDist (_pnts[_pnt_ids[1]], _pnts[_pnt_ids[2]]));
	if (tmp > _longest_edge) _longest_edge = tmp;
	tmp = MathLib::sqrDist (_pnts[_pnt_ids[0]], _pnts[_pnt_ids[2]]);
	if (tmp > _longest_edge) _longest_edge = tmp;
	_longest_edge = sqrt (_longest_edge);
}

void Triangle::setTriangle (size_t pnt_a, size_t pnt_b, size_t pnt_c)
{
	assert (pnt_a < _pnts.size() && pnt_b < _pnts.size() && pnt_c < _pnts.size());
	_pnt_ids[0] = pnt_a;
	_pnt_ids[1] = pnt_b;
	_pnt_ids[2] = pnt_c;

	_longest_edge = MathLib::sqrDist (_pnts[_pnt_ids[0]], _pnts[_pnt_ids[1]]);
	double tmp (MathLib::sqrDist (_pnts[_pnt_ids[1]], _pnts[_pnt_ids[2]]));
	if (tmp > _longest_edge) _longest_edge = tmp;
	tmp = MathLib::sqrDist (_pnts[_pnt_ids[0]], _pnts[_pnt_ids[2]]);
	if (tmp > _longest_edge) _longest_edge = tmp;
	_longest_edge = sqrt (_longest_edge);
}

bool Triangle::containsPoint (const double *pnt) const
{
	GEOLIB::Point const& a_tmp (*(_pnts[_pnt_ids[0]]));
	GEOLIB::Point const& b_tmp (*(_pnts[_pnt_ids[1]]));
	GEOLIB::Point const& c_tmp (*(_pnts[_pnt_ids[2]]));

	GEOLIB::Point s(a_tmp);
	for (size_t k(0); k<3; k++) {
		s[k] += b_tmp[k] + c_tmp[k];
		s[k] /= 3.0;
	}

	double eps (1e-2);
	GEOLIB::Point const a (a_tmp[0] + eps *(a_tmp[0]-s[0]),
			a_tmp[1] + eps *(a_tmp[1]-s[1]),
			a_tmp[2] + eps *(a_tmp[2]-s[2]));
	GEOLIB::Point const b (b_tmp[0] + eps *(b_tmp[0]-s[0]),
				b_tmp[1] + eps *(b_tmp[1]-s[1]),
				b_tmp[2] + eps *(b_tmp[2]-s[2]));
	GEOLIB::Point const c (c_tmp[0] + eps *(c_tmp[0]-s[0]),
				c_tmp[1] + eps *(c_tmp[1]-s[1]),
				c_tmp[2] + eps *(c_tmp[2]-s[2]));

//	GEOLIB::Point const& a (*(_pnts[_pnt_ids[0]]));
//	GEOLIB::Point const& b (*(_pnts[_pnt_ids[1]]));
//	GEOLIB::Point const& c (*(_pnts[_pnt_ids[2]]));

	const double delta (std::numeric_limits<double>::epsilon());
	const double upper (1+delta);

	// check special case where points of triangle have the same x-coordinate
	if (fabs(b[0]-a[0]) <= std::numeric_limits<double>::epsilon() &&
			fabs(c[0]-a[0]) <= std::numeric_limits<double>::epsilon()) {
		// all points of triangle have same x-coordinate
		if (fabs(pnt[0]-a[0]) / _longest_edge <= 1e-3) {
			// criterion: p-a = u0 * (b-a) + u1 * (c-a); 0 <= u0, u1 <= 1, u0+u1 <= 1
			MathLib::Matrix<double> mat (2,2);
			mat(0,0) = b[1] - a[1];
			mat(0,1) = c[1] - a[1];
			mat(1,0) = b[2] - a[2];
			mat(1,1) = c[2] - a[2];
			double y[2] = {pnt[1]-a[1], pnt[2]-a[2]};

			MathLib::GaussAlgorithm gauss (mat);
			gauss.execute (y);

			if (-delta <= y[0] && y[0] <= upper && -delta <= y[1] && y[1] <= upper
					&& y[0] + y[1] <= upper) {
				return true;
			} else {
				return false;
			}
		} else {
			return false;
		}
	}

	// check special case where points of triangle have the same y-coordinate
	if (fabs(b[1]-a[1]) <= std::numeric_limits<double>::epsilon() &&
			fabs(c[1]-a[1]) <= std::numeric_limits<double>::epsilon()) {
		// all points of triangle have same y-coordinate
		if (fabs(pnt[1]-a[1]) / _longest_edge <= 1e-3) {
			// criterion: p-a = u0 * (b-a) + u1 * (c-a); 0 <= u0, u1 <= 1, u0+u1 <= 1
			MathLib::Matrix<double> mat (2,2);
			mat(0,0) = b[0] - a[0];
			mat(0,1) = c[0] - a[0];
			mat(1,0) = b[2] - a[2];
			mat(1,1) = c[2] - a[2];
			double y[2] = {pnt[0]-a[0], pnt[2]-a[2]};

			MathLib::GaussAlgorithm gauss (mat);
			gauss.execute (y);

			if (-delta <= y[0] && y[0] <= upper && -delta <= y[1] && y[1] <= upper
					&& y[0] + y[1] <= upper) {
				return true;
			} else {
				return false;
			}
		} else {
			return false;
		}
	}

	// criterion: p-a = u0 * (b-a) + u1 * (c-a); 0 <= u0, u1 <= 1, u0+u1 <= 1
	MathLib::Matrix<double> mat (2,2);
	mat(0,0) = b[0] - a[0];
	mat(0,1) = c[0] - a[0];
	mat(1,0) = b[1] - a[1];
	mat(1,1) = c[1] - a[1];
	double y[2] = {pnt[0]-a[0], pnt[1]-a[1]};

	MathLib::GaussAlgorithm gauss (mat);
	gauss.execute (y);

	// check if the solution fulfills the third equation
	if (fabs((b[2]-a[2]) * y[0] + (c[2]-a[2]) * y[1] - (pnt[2] - a[2])) < 1e-3) {
		if (-delta <= y[0] && y[0] <= upper && -delta <= y[1] && y[1] <= upper &&
				y[0] + y[1] <= upper) {
			return true;
		}
		return false;
	} else {
		return false;
	}
}

} // end namespace GEOLIB
