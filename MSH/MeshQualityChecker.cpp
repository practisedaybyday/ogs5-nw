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
{
	if (_mesh) {
		_mesh_quality_messure.resize ((_mesh->getElementVector()).size(), 0);
	}
}

void MeshQualityChecker::check ()
{
	// get all elements of mesh
	const std::vector<Mesh_Group::CElem*>& msh_elem (_mesh->getElementVector());

	for (size_t k(0); k<msh_elem.size(); k++) {
		switch (msh_elem[k]->GetElementType()) {
		case MshElemType::TRIANGLE: {
			GEOLIB::Point* a (new GEOLIB::Point ((msh_elem[k]->GetNode(0))->getData ()));
			GEOLIB::Point* b (new GEOLIB::Point ((msh_elem[k]->GetNode(1))->getData ()));
			GEOLIB::Point* c (new GEOLIB::Point ((msh_elem[k]->GetNode(2))->getData ()));
			_mesh_quality_messure[k] = checkTriangle (a,b,c);
			delete a;
			delete b;
			delete c;
			break;
		}
		case MshElemType::QUAD: {
			GEOLIB::Point* a (new GEOLIB::Point ((msh_elem[k]->GetNode(0))->getData ()));
			GEOLIB::Point* b (new GEOLIB::Point ((msh_elem[k]->GetNode(1))->getData ()));
			GEOLIB::Point* c (new GEOLIB::Point ((msh_elem[k]->GetNode(2))->getData ()));
			GEOLIB::Point* d (new GEOLIB::Point ((msh_elem[k]->GetNode(3))->getData ()));
			_mesh_quality_messure[k] = checkQuad (a,b,c,d);
			delete a;
			delete b;
			delete c;
			delete d;
			break;
		}
		case MshElemType::TETRAHEDRON: {
			GEOLIB::Point* a (new GEOLIB::Point ((msh_elem[k]->GetNode(0))->getData ()));
			GEOLIB::Point* b (new GEOLIB::Point ((msh_elem[k]->GetNode(1))->getData ()));
			GEOLIB::Point* c (new GEOLIB::Point ((msh_elem[k]->GetNode(2))->getData ()));
			GEOLIB::Point* d (new GEOLIB::Point ((msh_elem[k]->GetNode(3))->getData ()));
			_mesh_quality_messure[k] = checkTetrahedron (a,b,c,d);
			delete a;
			delete b;
			delete c;
			delete d;
			break;
		}
//		case MshElemType::PRISM: {
//			GEOLIB::Point* a (new GEOLIB::Point ((msh_elem[k]->GetNode(0))->getData ()));
//			GEOLIB::Point* b (new GEOLIB::Point ((msh_elem[k]->GetNode(1))->getData ()));
//			GEOLIB::Point* c (new GEOLIB::Point ((msh_elem[k]->GetNode(2))->getData ()));
//			GEOLIB::Point* d (new GEOLIB::Point ((msh_elem[k]->GetNode(3))->getData ()));
//			GEOLIB::Point* e (new GEOLIB::Point ((msh_elem[k]->GetNode(4))->getData ()));
//			GEOLIB::Point* f (new GEOLIB::Point ((msh_elem[k]->GetNode(5))->getData ()));
//			_mesh_quality_messure[k] = checkPrism (a,b,c,d);
//			delete a;
//			delete b;
//			delete c;
//			delete d;
//			delete e;
//			delete f;
//			break;
//		}
		default:
			std::cout << "MeshQualityChecker::check () check for element type not implemented" << std::endl;
		}

	}
}

void MeshQualityChecker::getHistogramm (std::vector<size_t>& histogramm) const
{
	const size_t mesh_quality_messure_size (_mesh_quality_messure.size());
	const size_t histogramm_size (histogramm.size());
	for (size_t k(0); k<mesh_quality_messure_size; k++) {
		histogramm[static_cast<size_t>(_mesh_quality_messure[k] * histogramm_size)]++;
	}
}

double MeshQualityChecker::checkTriangle (GEOLIB::Point const * const a,
		GEOLIB::Point const * const b, GEOLIB::Point const * const c) const
{
	double len0 (sqrt(MATHLIB::sqrDist (b,a)));
	double len1 (sqrt(MATHLIB::sqrDist (b,c)));
	double len2 (sqrt(MATHLIB::sqrDist (a,c)));

	if (len0 < len1 && len0 < len2) {
		if (len1 < len2) {
			return len0/len2;
		} else {
			return len0/len1;
		}
	} else {
		if (len1 < len2) {
			if (len0 < len2) {
				return len1/len2;
			} else {
				return len1/len0;
			}
		} else {
			if (len0 < len1) {
				return len2/len1;
			} else {
				return len2/len0;
			}
		}
	}
}

double MeshQualityChecker::checkQuad (GEOLIB::Point const * const a, GEOLIB::Point const * const b,
		GEOLIB::Point const * const c, GEOLIB::Point const * const d) const
{
	double lengths[4] = {sqrt(MATHLIB::sqrDist (b,a)), sqrt(MATHLIB::sqrDist (c,b)),
			sqrt(MATHLIB::sqrDist (d,c)), sqrt(MATHLIB::sqrDist (a,d))};

	// sort lengths - since this is a very small array we use bubble sort
	for (size_t i(0); i<4; i++) {
		for (size_t j(i+1); j<4; j++) {
			if (lengths[i] >= lengths[j]) std::swap (lengths[i], lengths[j]);
		}
	}

	return lengths[0] / lengths[3];
}

double MeshQualityChecker::checkTetrahedron (GEOLIB::Point const * const a, GEOLIB::Point const * const b,
		GEOLIB::Point const * const c, GEOLIB::Point const * const d) const
{
	double lengths[6] = {sqrt(MATHLIB::sqrDist (b,a)), sqrt(MATHLIB::sqrDist (c,b)),
			sqrt(MATHLIB::sqrDist (c,a)), sqrt(MATHLIB::sqrDist (a,d)),
			sqrt(MATHLIB::sqrDist (b,d)), sqrt(MATHLIB::sqrDist (c,d))};

	// sort lengths - since this is a very small array we use bubble sort
	for (size_t i(0); i<6; i++) {
		for (size_t j(i+1); j<6; j++) {
			if (lengths[i] >= lengths[j]) std::swap (lengths[i], lengths[j]);
		}
	}

	return lengths[0] / lengths[5];
}

//double MeshQualityChecker::checkPrism (std::vector<GEOLIB::Point *> const & pnts) const
//{
//	double lengths[9] = {sqrt(MATHLIB::sqrDist (pnts[0],pnts[1])), sqrt(MATHLIB::sqrDist (pnts[1],pnts[2])),
//				sqrt(MATHLIB::sqrDist (pnts[2],pnts[0])), sqrt(MATHLIB::sqrDist (pnts[3],pnts[4])),
//				sqrt(MATHLIB::sqrDist (pnts[4],pnts[5])), sqrt(MATHLIB::sqrDist (pnts[5],pnts[3])),
//
//				sqrt(MATHLIB::sqrDist (pnts[4],pnts[5])), sqrt(MATHLIB::sqrDist (pnts[5],pnts[3])),
//				sqrt(MATHLIB::sqrDist (pnts[4],pnts[5]))};
//;
//}

}
