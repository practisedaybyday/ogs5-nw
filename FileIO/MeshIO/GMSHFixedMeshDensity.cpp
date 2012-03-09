/*
 * GMSHFixedMeshDensity.cpp
 *
 *  Created on: Mar 5, 2012
 *      Author: TF
 */

#include "MeshIO/GMSHFixedMeshDensity.h"

namespace FileIO {

GMSHFixedMeshDensity::GMSHFixedMeshDensity(double mesh_density) :
	_mesh_density(mesh_density)
{
}

void GMSHFixedMeshDensity::init(std::vector<GEOLIB::Point*> const& vec)
{
	// to avoid a warning here:
	(void)(vec);
}

void GMSHFixedMeshDensity::getMeshDensityAtPoint(GEOLIB::Point const*const pnt, std::ostream& out) const
{
	// to avoid a warning here:
	const_cast<GEOLIB::Point const*>(pnt);
	out << "," << _mesh_density;
}

} // end namespace FileIO
