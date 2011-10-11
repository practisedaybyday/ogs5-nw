/*
 * MeshQualityNormalisedVolumes.cpp
 *
 *  Created on: Mar 3, 2011
 *      Author: TF
 */

#include "MeshQualityNormalisedVolumes.h"

namespace MeshLib
{
MeshQualityNormalisedVolumes::MeshQualityNormalisedVolumes(
        CFEMesh const* const mesh) :
	MeshQualityChecker(mesh)
{}

void MeshQualityNormalisedVolumes::check()
{
	// get all elements of mesh
	const std::vector<MeshLib::CElem*>& msh_elem(_mesh->getElementVector());

	size_t error_count(0);
	double max_volume (0.0);
	double min_volume (std::numeric_limits<double>::max());

	for (size_t k(0); k < msh_elem.size(); k++)
	{
		MshElemType::type elem_type (msh_elem[k]->GetElementType());
		if (elem_type != MshElemType::LINE
		    && elem_type != MshElemType::TRIANGLE
		    && elem_type != MshElemType::QUAD)
		{
			double volume (msh_elem[k]->calcVolume());
			if (volume > max_volume)
				max_volume = volume;
			if (volume < sqrt(fabs(std::numeric_limits<double>::min())))
			{
				errorMsg(msh_elem[k], k);
				error_count++;
			}
			else if (volume < min_volume)
				min_volume = volume;
			_mesh_quality_messure[k] = volume;
		}
	}

	for (size_t k(0); k < msh_elem.size(); k++)
	{
		MshElemType::type elem_type (msh_elem[k]->GetElementType());
		if (elem_type != MshElemType::LINE
		    && elem_type != MshElemType::TRIANGLE
		    && elem_type != MshElemType::QUAD)
			_mesh_quality_messure[k] /= max_volume;
		else
			_mesh_quality_messure[k] = 1.1; // element has no valid value
	}

	std::cout << "MeshQualityNormalisedVolumes::check() min_volume: " << min_volume
	          << ", max_volume: " << max_volume << std::endl;
	if (error_count > 0)
		std::cout << "Warning: " << error_count << " elements with zero volume found." <<
		std::endl;
}

void MeshQualityNormalisedVolumes::getHistogramm (std::vector<size_t>& histogramm) const
{
	// get all elements of mesh
	std::vector<MeshLib::CElem*> const& msh_elem (_mesh->getElementVector());

	const size_t msh_elem_size (msh_elem.size());
	const size_t histogramm_size (histogramm.size() - 1);
	for (size_t k(0); k < msh_elem_size; k++)
	{
		MshElemType::type elem_type (msh_elem[k]->GetElementType());
		if (elem_type != MshElemType::LINE
		    && elem_type != MshElemType::TRIANGLE
		    && elem_type != MshElemType::QUAD)
			histogramm[static_cast<size_t>(_mesh_quality_messure[k] *
			                               histogramm_size)]++;
	}
}
} // end namespace MeshLib
