/*
 * LinearInterpolation.cpp
 *
 *  Created on: Sep 7, 2010
 *      Author: TF
 */

#include "LinearInterpolation.h"
#include "binarySearch.h"

#include <iostream>

namespace MATHLIB {

LinearInterpolation::LinearInterpolation(const std::vector<double>& supporting_points, const std::vector<double>& values_at_supp_pnts)
	: _supporting_points (supporting_points), _values_at_supp_pnts (values_at_supp_pnts)
{}

LinearInterpolation::LinearInterpolation(const std::vector<double>& supporting_points, const std::vector<double>& values_at_supp_pnts, const std::vector<double>& points_to_interpolate, std::vector<double>& values_at_interpol_pnts)
	: _supporting_points (supporting_points), _values_at_supp_pnts (values_at_supp_pnts)
{
	std::cout << "LinearInterpolation: input data:" << std::endl;
	for (size_t k(0); k<_supporting_points.size(); k++) {
		std::cout << "\t" << _supporting_points[k] << ", " << _values_at_supp_pnts[k] << std::endl;
	}
	values_at_interpol_pnts.clear();
	for (size_t k(0); k<points_to_interpolate.size(); k++)
		values_at_interpol_pnts.push_back (this->getValue (points_to_interpolate[k]));

	std::cout << "LinearInterpolation: results:" << std::endl;
	for (size_t k(0); k<points_to_interpolate.size(); k++) {
		std::cout << "\t" << points_to_interpolate[k] << ", " << values_at_interpol_pnts[k] << std::endl;
	}
	std::cout << "end LinearInterpolation" << std::endl;
}

LinearInterpolation::~LinearInterpolation()
{}

double LinearInterpolation::getValue ( double pnt_to_interpolate )
{
	if (pnt_to_interpolate < _supporting_points[0] || _supporting_points[_supporting_points.size()-1] < pnt_to_interpolate)
		return std::numeric_limits<double>::min ();

	// search interval that has the point inside
	size_t interval_idx (std::numeric_limits<size_t>::max());
	for (size_t k(1); k<_supporting_points.size() && interval_idx == std::numeric_limits<size_t>::max(); k++) {
		if (_supporting_points[k-1] <= pnt_to_interpolate && pnt_to_interpolate <= _supporting_points[k]) {
			interval_idx = k-1;
		}
	}

	// compute linear interpolation polynom: y = m * x + n
	double m ((_values_at_supp_pnts[interval_idx+1] - _values_at_supp_pnts[interval_idx]) / (_supporting_points[interval_idx+1] - _supporting_points[interval_idx]));
	double n (_values_at_supp_pnts[interval_idx+1] - m * _supporting_points[interval_idx+1]);

	return m * pnt_to_interpolate + n;
}

} // end MATHLIB
