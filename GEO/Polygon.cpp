/*
 * Polygon.cpp
 *
 *  Created on: Jun 21, 2010
 *      Author: TF
 */

#include "Polygon.h"

// MathLib
#include "AnalyticalGeometry.h"
#include "MathTools.h"

// Base
#include "quicksort.h"

namespace GEOLIB {

Polygon::Polygon(const Polyline &ply) :
	Polyline (ply),
	_maxx (std::numeric_limits<double>::min()),
	_maxy (std::numeric_limits<double>::min())
{
	size_t n_nodes (getNumberOfPoints()-1);
	for (size_t k(0); k<n_nodes; k++) {
		if ((*(getPoint(k)))[0] > _maxx) _maxx = (*(getPoint(k)))[0];
		if ((*(getPoint(k)))[1] > _maxy) _maxy = (*(getPoint(k)))[1];
	}

	double s0 (0.1), s1(0.2); // small random values
	_maxx += _maxx*s0;
	_maxy += _maxy*s1;
}

Polygon::~Polygon()
{}

bool Polygon::isPntInPolygon (const GEOLIB::Point& pnt) const
{
	size_t n_intersections (0);
	GEOLIB::Point s;

	if (_simple_polygon_list.empty ()) {
		const size_t n_nodes (getNumberOfPoints()-1);
		const GEOLIB::Point other (_maxx, _maxy, pnt[2]);
		for (size_t k(0); k<n_nodes; k++) {
			if (MATHLIB::lineSegmentIntersect (*(getPoint(k)), *(getPoint(k+1)), other, pnt, s)) {
				n_intersections++;
			}
		}
	} else {
		for (std::list<Polygon*>::const_iterator it (_simple_polygon_list.begin());
			it != _simple_polygon_list.end(); ++it) {
			const Polygon* polygon (*it);
			const GEOLIB::Point other (_maxx, _maxy, pnt[2]);
			const size_t n_nodes_simple_polygon (polygon->getNumberOfPoints()-1);
			for (size_t k(0); k<n_nodes_simple_polygon; k++) {
				if (MATHLIB::lineSegmentIntersect (*(polygon->getPoint(k)), *(polygon->getPoint(k+1)), other, pnt, s)) {
					n_intersections++;
				}
			}
		}
	}
	if (n_intersections%2 == 1) return true;

	// check if point is at the border
	if (_simple_polygon_list.empty ()) {
		const size_t n_nodes (getNumberOfPoints());
		for (size_t k(0); k<n_nodes; k++) {
			if (MATHLIB::sqrDist (getPoint(k), &pnt) < std::numeric_limits<double>::min())
				return true;
		}
	} else {
		for (std::list<Polygon*>::const_iterator it (_simple_polygon_list.begin());
			it != _simple_polygon_list.end(); ++it) {
			const size_t n_nodes_simple_polygon ((*it)->getNumberOfPoints());
			for (size_t k(0); k<n_nodes_simple_polygon; k++) {
				if (MATHLIB::sqrDist ((*it)->getPoint(k), &pnt)) {
					return true;
				}
			}
		}
	}

	return false;
}

bool Polygon::isPolylineInPolygon (const Polyline& ply) const
{
	size_t ply_size (ply.getNumberOfPoints()), cnt (0);
	for (size_t k(0); k<ply_size; k++) {
		if (isPntInPolygon (*(ply[k]))) {
			cnt++;
		}
	}
	if (cnt == ply_size)
		return true;
	return false;
}

const std::list<Polygon*>& Polygon::getListOfSimplePolygons()
{
	if (_simple_polygon_list.empty())
		_simple_polygon_list.push_back (this);
	return _simple_polygon_list;
}

void Polygon::computeListOfSimplePolygons ()
{
	_simple_polygon_list.push_back (this);
	splitPolygonAtPoint (_simple_polygon_list.begin());
	splitPolygonAtIntersection (_simple_polygon_list.begin());
}

void Polygon::splitPolygonAtIntersection (std::list<Polygon*>::iterator polygon_it)
{
	size_t idx0 (0), idx1 (0);
	while (polygon_it != _simple_polygon_list.end()) {
		GEOLIB::Point *intersection_pnt (new GEOLIB::Point);
		bool is_simple (!MATHLIB::lineSegmentsIntersect (*polygon_it, idx0, idx1, *intersection_pnt));
		if (!is_simple) {
			// adding intersection point to pnt_vec
			size_t intersection_pnt_id (_ply_pnts.size());
			const_cast<std::vector<Point*>& >(_ply_pnts).push_back (intersection_pnt);

			// split Polygon
			if (idx0 > idx1) BASELIB::swap (idx0, idx1);

			GEOLIB::Polygon* polygon0 (new GEOLIB::Polygon((*polygon_it)->getPointsVec()));
			for (size_t k(0); k<=idx0; k++) polygon0->addPoint ((*polygon_it)->getPointID (k));
			polygon0->addPoint (intersection_pnt_id);
			for (size_t k(idx1+1); k<(*polygon_it)->getNumberOfPoints(); k++)
				polygon0->addPoint ((*polygon_it)->getPointID (k));

			GEOLIB::Polygon* polygon1 (new GEOLIB::Polygon((*polygon_it)->getPointsVec()));
			polygon1->addPoint (intersection_pnt_id);
			for (size_t k(idx0+1); k<=idx1; k++) polygon1->addPoint ((*polygon_it)->getPointID (k));
			polygon1->addPoint (intersection_pnt_id);

			// remove original polyline and add two new polylines
			std::list<GEOLIB::Polygon*>::iterator polygon0_it, polygon1_it;
			polygon_it = _simple_polygon_list.erase (polygon_it);
			polygon1_it = _simple_polygon_list.insert (polygon_it, polygon1);
			polygon0_it = _simple_polygon_list.insert (polygon1_it, polygon0);

			splitPolygonAtIntersection (polygon0_it);
			splitPolygonAtIntersection (polygon1_it);
		} else {
			delete intersection_pnt;
		}
		++polygon_it;
	}
}

void Polygon::splitPolygonAtPoint (std::list<GEOLIB::Polygon*>::iterator polygon_it)
{
	size_t n ((*polygon_it)->getNumberOfPoints()-1), idx0 (0), idx1(0);
	size_t *id_vec (new size_t[n]), *perm (new size_t[n]);
	for (size_t k(0); k<n; k++) {
		id_vec[k] = (*polygon_it)->getPointID (k);
		perm[k] = k;
	}

	quicksort (id_vec, 0, n, perm);

	for (size_t k(0); k<n-1; k++) {
		if (id_vec[k] == id_vec[k+1]) {
			idx0 = perm[k];
			idx1 = perm[k+1];
			delete [] perm;
			delete [] id_vec;

			if (idx0 > idx1) BASELIB::swap (idx0, idx1);

			GEOLIB::Polygon* polygon0 (new GEOLIB::Polygon((*polygon_it)->getPointsVec()));
			for (size_t k(0); k<=idx0; k++) polygon0->addPoint ((*polygon_it)->getPointID (k));
			for (size_t k(idx1+1); k<(*polygon_it)->getNumberOfPoints(); k++) polygon0->addPoint ((*polygon_it)->getPointID (k));

			GEOLIB::Polygon* polygon1 (new GEOLIB::Polygon((*polygon_it)->getPointsVec()));
			for (size_t k(idx0); k<=idx1; k++) polygon1->addPoint ((*polygon_it)->getPointID (k));

			// remove original polygon and add two new polygons
			std::list<GEOLIB::Polygon*>::iterator polygon0_it, polygon1_it;
			polygon1_it = _simple_polygon_list.insert (_simple_polygon_list.erase (polygon_it), polygon1);
			polygon0_it = _simple_polygon_list.insert (polygon1_it, polygon0);

			splitPolygonAtPoint (polygon0_it);
			splitPolygonAtPoint (polygon1_it);

			return;
		}
	}
	delete [] perm;
	delete [] id_vec;
}

} // end namespace GEOLIB
