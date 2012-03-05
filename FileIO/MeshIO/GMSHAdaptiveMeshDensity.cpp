/*
 * GMSHAdaptiveMeshDensity.cpp
 *
 *  Created on: Mar 5, 2012
 *      Author: TF
 */

#include "MeshIO/GMSHAdaptiveMeshDensity.h"

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
	_quad_tree = new GEOLIB::QuadTree<GEOLIB::Point> (min, max, _max_pnt_per_leaf);
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

std::ostream& GMSHAdaptiveMeshDensity::getMeshDensityAtPoint(GEOLIB::Point const*const, std::ostream&)
{
	;
}


} // end namespace FileIO
