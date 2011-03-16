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
			std::cout << "Error in MeshQualityNormalisedVolumes::check() - Volume of element is below double precision minimum value." << std::endl;
			std::cout << "Points of " << MshElemType2String(msh_elem[k]->GetElementType()) << "-Element " << k << ": " << std::endl;
			for (int i(0); i<msh_elem[k]->GetVertexNumber(); i++)
				std::cout << "\t Node " << i << " " << GEOLIB::Point((msh_elem[k]->GetNode(i))->getData()) << std::endl;
		}
		_mesh_quality_messure[k] = volume;
	}
	for (size_t k(0); k < msh_elem.size(); k++) {
		_mesh_quality_messure[k] /= max_volume;
	}
}

} // end namespace Mesh_Group
