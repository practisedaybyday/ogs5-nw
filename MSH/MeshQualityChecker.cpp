/*
 * MeshQualityChecker.cpp
 *
 *  Created on: Dec 8, 2010
 *      Author: TF
 */

#include "MeshQualityChecker.h"
#include "msh_elem.h"

namespace Mesh_Group {

MeshQualityChecker::MeshQualityChecker(CFEMesh const * const mesh) :
	_mesh (mesh), _static_histogramm (100, 0)
{}

void MeshQualityChecker::checkTriangles ()
{
	// get all elements of mesh
	const std::vector<Mesh_Group::CElem*>& msh_elem (_mesh->getElementVector());

	double limit (0.26);

	for (size_t k(0); k<msh_elem.size(); k++) {
		if (msh_elem[k]->GetElementType() == MshElemType::TRIANGLE) {
			const double* a ((msh_elem[k]->GetNode(0))->getData ());
			const double* b ((msh_elem[k]->GetNode(1))->getData ());
			const double* c ((msh_elem[k]->GetNode(2))->getData ());

			double ab (sqrt(MATHLIB::sqrDist (b,a)));
			double bc (sqrt(MATHLIB::sqrDist (b,c)));
			double ca (sqrt(MATHLIB::sqrDist (a,c)));

			if (ab < bc && ab < ca) {
				if (bc < ca) {
					_static_histogramm[static_cast<size_t>(ab/ca * 100)]++;
					if (ab/ca < limit) std::cout << k << std::endl;
				} else {
					_static_histogramm[static_cast<size_t>(ab/bc * 100)]++;
					if (ab/bc < limit) std::cout << k << std::endl;
				}
			} else {
				if (bc < ca) {
					if (ab < ca) {
						_static_histogramm[static_cast<size_t>(bc/ca * 100)]++;
						if (bc/ca < limit) std::cout << k << std::endl;
					} else {
						_static_histogramm[static_cast<size_t>(bc/ab * 100)]++;
						if (bc/ab < limit) std::cout << k << std::endl;
					}
				} else {
					if (ab < bc) {
						_static_histogramm[static_cast<size_t>(ca/bc * 100)]++;
						if (ca/bc < limit) std::cout << k << std::endl;
					} else {
						_static_histogramm[static_cast<size_t>(ca/ab * 100)]++;
						if (ca/ab < limit) std::cout << k << std::endl;
					}
				}
			}
		}
	}
}

MeshQualityChecker::~MeshQualityChecker()
{}

}
