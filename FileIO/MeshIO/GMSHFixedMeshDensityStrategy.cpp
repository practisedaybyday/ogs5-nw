/*
 * GMSHFixedMeshDensityStrategy.cpp
 *
 *  Created on: Mar 5, 2012
 *      Author: TF
 */

#include "MeshIO/GMSHFixedMeshDensityStrategy.h"

namespace FileIO {

GMSHFixedMeshDensityStrategy::GMSHFixedMeshDensityStrategy(double mesh_density) :
	_mesh_density(mesh_density)
{
}

void GMSHFixedMeshDensityStrategy::init(std::vector<GEOLIB::Point*> const& vec)
{
	// to avoid a warning here:
	(void*)(vec);
}

std::ostream& GMSHFixedMeshDensityStrategy::getMeshDensityAtPoint(GEOLIB::Point const*const pnt, std::ostream& out) const
{
	// to avoid a warning here:
	const_cast<GEOLIB::Point const*>(pnt);

	out << _mesh_density;
	return out;
}

} // end namespace FileIO
