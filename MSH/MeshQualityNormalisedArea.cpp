/*
 * MeshQualityNormalisedArea.cpp
 *
 * 2011/03/17 KR Initial Implementation
 */

#include "MeshQualityNormalisedArea.h"
#include "mathlib.cpp"

namespace Mesh_Group {

MeshQualityNormalisedArea::MeshQualityNormalisedArea(CFEMesh const * const mesh)
: MeshQualityChecker(mesh)
{}

void MeshQualityNormalisedArea::check()
{
	// get all elements of mesh
	const std::vector<Mesh_Group::CElem*>& msh_elem(_mesh->getElementVector());

	double max_volume (0.0);

	size_t nElems(msh_elem.size());
	for (size_t k(0); k < nElems; k++)
	{
		MshElemType::type elem_type (msh_elem[k]->GetElementType());
		if (elem_type != MshElemType::LINE)
		{
			double min_volume(std::numeric_limits<double>::max());
			if (elem_type == MshElemType::TRIANGLE || elem_type != MshElemType::QUAD)
			{
				min_volume = msh_elem[k]->calcVolume();
				if (min_volume < sqrt(fabs(std::numeric_limits<double>::min()))) errorMsg(msh_elem[k], k);
			}
			else
			{
				size_t nFaces (msh_elem[k]->GetFacesNumber());

				int face_node_index[4];
				for (size_t i=0; i<nFaces; i++)
				{
					size_t nNodes = msh_elem[k]->GetElementFaceNodes(i, face_node_index);

					double volume(0);
					if (nNodes==3)		// face is a triangle
						volume = ComputeDetTri(msh_elem[k]->GetNode(0)->getData(), msh_elem[k]->GetNode(1)->getData(), msh_elem[k]->GetNode(2)->getData());
					else if (nNodes==4)	// face is a quad
						volume = ComputeDetTri(msh_elem[k]->GetNode(0)->getData(), msh_elem[k]->GetNode(1)->getData(), msh_elem[k]->GetNode(2)->getData())
							   + ComputeDetTri(msh_elem[k]->GetNode(2)->getData(), msh_elem[k]->GetNode(3)->getData(), msh_elem[k]->GetNode(0)->getData());
					else std::cout << "Error in MeshQualityNormalisedArea::check()" << std::endl;

					if (volume < sqrt(fabs(std::numeric_limits<double>::min()))) errorMsg(msh_elem[k], k);
					if (volume < min_volume) min_volume = volume;
				}
			}
			if (min_volume > max_volume) max_volume = min_volume;
			_mesh_quality_messure[k] = min_volume;
		}
	}

	for (size_t k(0); k < nElems; k++) {
		MshElemType::type elem_type (msh_elem[k]->GetElementType());
		_mesh_quality_messure[k] = (elem_type == MshElemType::LINE) ? 1.1 : (_mesh_quality_messure[k] / max_volume);
	}
}

} // end namespace Mesh_Group
