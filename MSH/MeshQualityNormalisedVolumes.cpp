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
	MeshQualityVolumes(mesh)
{}

void MeshQualityNormalisedVolumes::check()
{
	MeshQualityVolumes::check();

	// Normalisation.
	typedef std::vector<QualityType>::iterator QI;
	for (QI q = _mesh_quality_measure.begin();
		q != _mesh_quality_measure.end(); ++q)
	{
		if (*q)
			**q /= _maximum;
	}
}

} // end namespace MeshLib
