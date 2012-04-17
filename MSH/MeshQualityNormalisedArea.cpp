/*
 * MeshQualityNormalisedArea.cpp
 *
 * 2011/03/17 KR Initial Implementation
 */

#include "MeshQualityNormalisedArea.h"

namespace MeshLib
{
MeshQualityNormalisedArea::MeshQualityNormalisedArea(CFEMesh const* const mesh)
	: MeshQualityArea(mesh)
{}

void MeshQualityNormalisedArea::check()
{
	MeshQualityArea::check();

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
