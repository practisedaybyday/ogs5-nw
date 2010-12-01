/*
 * Polyline.cpp
 *
 *  Created on: Jun 21, 2010
 *      Author: TF
 */

// Base
#include "swap.h"

#include "Polyline.h"

namespace GEOLIB {

Polyline::Polyline(const std::vector<Point*>& pnt_vec) :
	GeoObject(), _ply_pnts(pnt_vec)
{}

Polyline::Polyline(const Polyline& ply) :
	GeoObject(), _ply_pnts (ply._ply_pnts)
{
	for (size_t k(0); k<ply.getNumberOfPoints(); ++k) {
		_ply_pnt_ids.push_back (ply.getPointID (k));
	}

	if (ply.getNumberOfPoints() > 0) {
		for (size_t k(0); k<ply.getNumberOfPoints()-1; ++k) {
			_length.push_back (ply.getLength (k));
		}
	}
}

void Polyline::write(std::ostream &os) const
{
	size_t size(_ply_pnt_ids.size());
	for (size_t k(0); k < size; k++) {
		os << *(_ply_pnts[_ply_pnt_ids[k]]) << std::endl;
	}
}

void Polyline::addPoint(size_t point_id)
{
	assert(point_id < _ply_pnts.size());
	size_t n_pnts (_ply_pnt_ids.size());
	_ply_pnt_ids.push_back(point_id);

	if (n_pnts > 0) {
		double act_dist (sqrt(MATHLIB::sqrDist (_ply_pnts[_ply_pnt_ids[n_pnts-1]], _ply_pnts[point_id])));
		double dist_until_now (0.0);
		if (n_pnts > 1)
			dist_until_now = _length[n_pnts - 2];

		_length.push_back (dist_until_now + act_dist);
	}
}

size_t Polyline::getNumberOfPoints() const
{
	return _ply_pnt_ids.size();
}

bool Polyline::isClosed() const
{
	if (_ply_pnt_ids.front() == _ply_pnt_ids.back())
		return true;
	else
		return false;
}

size_t Polyline::getPointID(size_t i) const
{
	assert(i < _ply_pnt_ids.size());
	return _ply_pnt_ids[i];
}

const Point* Polyline::operator[](size_t i) const
{
	assert(i < _ply_pnt_ids.size());
	return _ply_pnts[_ply_pnt_ids[i]];
}

const Point* Polyline::getPoint(size_t i) const
{
	assert(i < _ply_pnt_ids.size());
 	return _ply_pnts[_ply_pnt_ids[i]];
}

const std::vector<Point*> & Polyline::getPointsVec () const
{
	return _ply_pnts;
}

double Polyline::getLength (size_t k) const
{
	assert(k < _length.size());
	return _length[k];
}

const std::vector<double>& Polyline::getLengthVec () const
{
	return _length;
}

std::ostream& operator<< (std::ostream &os, const Polyline &pl)
{
	pl.write (os);
	return os;
}

bool containsEdge (const Polyline& ply, size_t id0, size_t id1)
{
	if (id0 == id1) {
		std::cerr << "no valid edge id0 == id1 == " << id0 << std::endl;
		return false;
	}
	if (id0 > id1) BASELIB::swap (id0,id1);
	const size_t n (ply.getNumberOfPoints() - 1);
	for (size_t k(0); k<n; k++) {
		size_t ply_pnt0 (ply.getPointID (k));
		size_t ply_pnt1 (ply.getPointID (k+1));
		if (ply_pnt0 > ply_pnt1)
			BASELIB::swap (ply_pnt0, ply_pnt1);
		if (ply_pnt0 == id0 && ply_pnt1 == id1)
			return true;
	}
	return false;
}

} // end namespace GEOLIB
