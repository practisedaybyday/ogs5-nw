/*
 * MeshQualityChecker.cpp
 *
 *  Created on: Dec 8, 2010
 *      Author: TF
 */

#include "MeshQualityChecker.h"
#include "msh_elem.h"
#include <cmath>

namespace MeshLib
{
MeshQualityChecker::MeshQualityChecker(CFEMesh const* const mesh) :
	_mesh (mesh), _static_histogramm (100, 0)
{
	if (_mesh)
		_mesh_quality_messure.resize ((_mesh->getElementVector()).size(), 0);
}

void MeshQualityChecker::getHistogramm (std::vector<size_t>& histogramm) const
{
	// get all elements of mesh
	const std::vector<MeshLib::CElem*>& msh_elem (_mesh->getElementVector());

	const size_t msh_elem_size (msh_elem.size());
	const size_t histogramm_size (histogramm.size() - 1);
	for (size_t k(0); k < msh_elem_size; k++)
		if (msh_elem[k]->GetElementType() != MshElemType::LINE)
			histogramm[static_cast<size_t>(_mesh_quality_messure[k] *
			                               histogramm_size)]++;

}

void MeshQualityChecker::errorMsg (CElem* elem, size_t idx) const
{
	std::cout <<
	"Error in MeshQualityChecker::check() - Calculated value of element is below double precision minimum."
	          << std::endl;
	std::cout << "Points of " << MshElemType2String(elem->GetElementType()) << "-Element " <<
	idx << ": " << std::endl;
	for (int i(0); i < elem->GetVertexNumber(); i++)
		std::cout << "\t Node " << i << " " <<
		GEOLIB::Point((elem->GetNode(i))->getData()) << std::endl;
}
}
