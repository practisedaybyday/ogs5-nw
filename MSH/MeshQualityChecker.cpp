/*
 * MeshQualityChecker.cpp
 *
 *  Created on: Dec 8, 2010
 *      Author: TF
 */

#include "MeshQualityChecker.h"
#include "msh_elem.h"
#include <cmath>

namespace Mesh_Group {

MeshQualityChecker::MeshQualityChecker(CFEMesh const * const mesh) :
	_mesh (mesh), _static_histogramm (100, 0)
{
	if (_mesh) {
		_mesh_quality_messure.resize ((_mesh->getElementVector()).size(), 0);
	}
}

void MeshQualityChecker::getHistogramm (std::vector<size_t>& histogramm) const
{
	// get all elements of mesh
	const std::vector<Mesh_Group::CElem*>& msh_elem (_mesh->getElementVector());

	const size_t msh_elem_size (msh_elem.size());
	const size_t histogramm_size (histogramm.size()-1);
	for (size_t k(0); k<msh_elem_size; k++) {
		if (msh_elem[k]->GetElementType() != MshElemType::LINE) {
			histogramm[static_cast<size_t>(_mesh_quality_messure[k] * histogramm_size)]++;
		}
	}
}

}
