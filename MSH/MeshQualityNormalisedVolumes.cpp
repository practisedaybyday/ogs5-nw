/*
 * MeshQualityNormalisedVolumes.cpp
 *
 *  Created on: Mar 3, 2011
 *      Author: TF
 */

#include "MeshQualityNormalisedVolumes.h"

namespace Mesh_Group {

MeshQualityNormalisedVolumes::MeshQualityNormalisedVolumes(
		CFEMesh const * const mesh) :
	MeshQualityChecker(mesh)
{}

void MeshQualityNormalisedVolumes::check()
{
	// get all elements of mesh
	const std::vector<Mesh_Group::CElem*>& msh_elem(_mesh->getElementVector());

	double max_volume (0.0);

	for (size_t k(0); k < msh_elem.size(); k++) {
		double volume (msh_elem[k]->calcVolume());
		if (volume > max_volume) max_volume = volume;
		if (volume < sqrt(fabs(std::numeric_limits<double>::min()))) {
			std::cout << "volume of element " << k << " is " << volume << std::endl;
			std::cout << "points of tetrahedron: " << std::endl;
			std::cout << "\t" << 0 << " " << GEOLIB::Point((msh_elem[k]->GetNode(0))->getData()) << std::endl;
			std::cout << "\t" << 1 << " " << GEOLIB::Point((msh_elem[k]->GetNode(1))->getData()) << std::endl;
			std::cout << "\t" << 2 << " " << GEOLIB::Point((msh_elem[k]->GetNode(2))->getData()) << std::endl;
			std::cout << "\t" << 3 << " " << GEOLIB::Point((msh_elem[k]->GetNode(3))->getData()) << std::endl;
		}
		_mesh_quality_messure[k] = volume;
	}
	for (size_t k(0); k < msh_elem.size(); k++) {
		_mesh_quality_messure[k] /= max_volume;
	}
}

} // end namespace Mesh_Group
