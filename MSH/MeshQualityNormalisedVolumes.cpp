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
		MshElemType::type elem_type (msh_elem[k]->GetElementType());
		if (elem_type != MshElemType::LINE
				&& elem_type != MshElemType::TRIANGLE
				&& elem_type != MshElemType::QUAD) {
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
	}

	for (size_t k(0); k < msh_elem.size(); k++) {
		MshElemType::type elem_type (msh_elem[k]->GetElementType());
		if (elem_type != MshElemType::LINE
			&& elem_type != MshElemType::TRIANGLE
			&& elem_type != MshElemType::QUAD) {
			_mesh_quality_messure[k] /= max_volume;
		} else {
			_mesh_quality_messure[k] = 1.1; // element has no valid value
		}
	}
}

void MeshQualityNormalisedVolumes::getHistogramm (std::vector<size_t>& histogramm) const
{
	// get all elements of mesh
	std::vector<Mesh_Group::CElem*> const& msh_elem (_mesh->getElementVector());

	const size_t msh_elem_size (msh_elem.size());
	const size_t histogramm_size (histogramm.size()-1);
	for (size_t k(0); k<msh_elem_size; k++) {
		MshElemType::type elem_type (msh_elem[k]->GetElementType());
		if (elem_type != MshElemType::LINE
				&& elem_type != MshElemType::TRIANGLE
				&& elem_type != MshElemType::QUAD) {
			histogramm[static_cast<size_t>(_mesh_quality_messure[k] * histogramm_size)]++;
		}
	}
}

} // end namespace Mesh_Group
