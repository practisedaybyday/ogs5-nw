/*
 * GMSHAdaptiveMeshDensity.cpp
 *
 *  Created on: Mar 5, 2012
 *      Author: TF
 */

#include <list>

// FileIO
#include "MeshIO/GMSHAdaptiveMeshDensity.h"

// GEOLIB
#include "Polygon.h"

namespace FileIO {

GMSHAdaptiveMeshDensity::GMSHAdaptiveMeshDensity(double pnt_density, double station_density,
				size_t max_pnts_per_leaf) :
	_pnt_density(pnt_density), _station_density(station_density),
	_max_pnts_per_leaf(max_pnts_per_leaf), _quad_tree(NULL)
{
}

GMSHAdaptiveMeshDensity::~GMSHAdaptiveMeshDensity()
{
	delete _quad_tree;
}

void GMSHAdaptiveMeshDensity::init(std::vector<GEOLIB::Point*> const& pnts)
{
	// *** QuadTree - determining bounding box
#ifndef NDEBUG
	std::cout << "[GMSHAdaptiveMeshDensity::init]" << std::endl;
	std::cout << "\tcomputing axis aligned bounding box (2D) for quadtree ... " << std::flush;
#endif
	GEOLIB::Point min(pnts[0]->getData()), max(pnts[0]->getData());
	size_t n_pnts(pnts.size());
	for (size_t k(1); k<n_pnts; k++) {
		for (size_t j(0); j<2; j++)
			if ((*(pnts[k]))[j] < min[j]) min[j] = (*(pnts[k]))[j];
		for (size_t j(0); j<2; j++)
			if ((*(pnts[k]))[j] > max[j]) max[j] = (*(pnts[k]))[j];
	}
	min[2] = 0.0;
	max[2] = 0.0;
#ifndef NDEBUG
	std::cout << "ok" << std::endl;
#endif

	// *** QuadTree - create object
#ifndef NDEBUG
	std::cout << "\tcreating quadtree ... " << std::flush;
#endif
	_quad_tree = new GEOLIB::QuadTree<GEOLIB::Point> (min, max, _max_pnts_per_leaf);
#ifndef NDEBUG
	std::cout << "ok" << std::endl;
#endif

	// *** QuadTree - insert points
	addPoints(pnts);
}

void GMSHAdaptiveMeshDensity::addPoints(std::vector<GEOLIB::Point*> const& pnts)
{
	// *** QuadTree - insert points
	const size_t n_pnts(pnts.size());
#ifndef NDEBUG
	std::cout << "\tinserting " << n_pnts << " points into quadtree ... " <<
	std::flush;
#endif
	for (size_t k(0); k < n_pnts; k++)
		_quad_tree->addPoint(pnts[k]);
#ifndef NDEBUG
	std::cout << "ok" << std::endl;
#endif
}

void GMSHAdaptiveMeshDensity::getMeshDensityAtPoint(GEOLIB::Point const*const pnt, std::ostream& out) const
{
	GEOLIB::Point ll, ur;
	_quad_tree->getLeaf(*pnt, ll, ur);
	out << "," << (_pnt_density * (ur[0] - ll[0]));
}

void GMSHAdaptiveMeshDensity::getMeshDensityAtStation(GEOLIB::Point const*const pnt, std::ostream& out)const
{
	GEOLIB::Point ll, ur;
	_quad_tree->getLeaf(*pnt, ll, ur);
	out << "," << (_station_density * (ur[0] - ll[0]));
}

void GMSHAdaptiveMeshDensity::writeSteinerPoints(GEOLIB::Polygon const*const bounding_polygon,
				size_t &pnt_idx_offset, size_t sfc_idx_offset, std::ostream & out) const
{
	// write Steiner points
	std::list<GEOLIB::QuadTree<GEOLIB::Point>*> leaf_list;
	_quad_tree->getLeafs(leaf_list);

	for (std::list<GEOLIB::QuadTree<GEOLIB::Point>*>::const_iterator it(leaf_list.begin()); it
					!= leaf_list.end(); it++) {
		if ((*it)->getPoints().empty()) {
			// compute point from square
			GEOLIB::Point ll, rr;
			(*it)->getSquarePoints(ll, rr);
			GEOLIB::Point mid_point(0.5 * (rr[0] + ll[0]), 0.5 * (rr[1] + ll[1]), 0.5 * (rr[2]
							+ ll[2]));
			if (bounding_polygon->isPntInPolygon(mid_point)) {
				out << "Point(" << pnt_idx_offset << ") = {" << mid_point[0] << "," << mid_point[1]
								<< "," << mid_point[2] << "," << 0.5 * (rr[0] - ll[0]) << "};"
								<< std::endl;
				out << "Point {" << pnt_idx_offset << "} In Surface {" << sfc_idx_offset - 1 << "};"
								<< std::endl;
				pnt_idx_offset++;

			}
		}
	}
}

} // end namespace FileIO
