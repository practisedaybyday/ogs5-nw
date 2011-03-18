/*
 * MeshQualityEquiAngleSkew.cpp
 *
 *  Created on: Mar 17, 2011
 *      Author: TF
 */

#include "MeshQualityEquiAngleSkew.h"

// MathLib
#include "MathTools.h"

#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace Mesh_Group {

MeshQualityEquiAngleSkew::MeshQualityEquiAngleSkew(CFEMesh const * const mesh) :
	MeshQualityChecker(mesh), M_PI_THIRD (M_PI/3.0), M_PI_HALF (M_PI/2.0)
{}

MeshQualityEquiAngleSkew::~MeshQualityEquiAngleSkew()
{}

void MeshQualityEquiAngleSkew::check ()
{
	// get all elements of mesh
	const std::vector<Mesh_Group::CElem*>& msh_elem(_mesh->getElementVector());

	for (size_t k(0); k < msh_elem.size(); k++) {
		switch (msh_elem[k]->GetElementType()) {
		case MshElemType::TRIANGLE:
			_mesh_quality_messure[k] = checkTriangle (msh_elem[k]);
			break;
		case MshElemType::QUAD:
			_mesh_quality_messure[k] = checkQuad (msh_elem[k]);
			break;
		case MshElemType::TETRAHEDRON:
			_mesh_quality_messure[k] = checkTetrahedron (msh_elem[k]);
			break;
		case MshElemType::HEXAHEDRON:
			_mesh_quality_messure[k] = checkHexahedron (msh_elem[k]);
			break;
		}

	}
}

double MeshQualityEquiAngleSkew::checkTriangle (CElem const * const elem) const
{
	double const * const node0 (elem->GetNode(0)->getData());
	double const * const node1 (elem->GetNode(1)->getData());
	double const * const node2 (elem->GetNode(2)->getData());

	double min_angle (MATHLIB::getAngle (node2, node0, node1));
	double max_angle (min_angle);

	double angle (MATHLIB::getAngle (node0, node1, node2));
	if (angle < min_angle) min_angle = angle;
	if (angle > max_angle) max_angle = angle;

	angle = MATHLIB::getAngle (node1, node2, node0);
	if (angle < min_angle) min_angle = angle;
	if (angle > max_angle) max_angle = angle;

	return 1.0 - std::max((max_angle - M_PI_THIRD)/(M_PI - M_PI_THIRD), (M_PI_THIRD - min_angle)/(M_PI_THIRD));
}

double MeshQualityEquiAngleSkew::checkQuad (CElem const * const elem) const
{
	double const * const node0 (elem->GetNode(0)->getData());
	double const * const node1 (elem->GetNode(1)->getData());
	double const * const node2 (elem->GetNode(2)->getData());
	double const * const node3 (elem->GetNode(3)->getData());

	double min_angle (MATHLIB::getAngle (node3, node0, node1));
	double max_angle (min_angle);

	double angle (MATHLIB::getAngle (node0, node1, node2));
	if (angle < min_angle) min_angle = angle;
	if (angle > max_angle) max_angle = angle;

	angle = MATHLIB::getAngle (node1, node2, node3);
	if (angle < min_angle) min_angle = angle;
	if (angle > max_angle) max_angle = angle;

	angle = MATHLIB::getAngle (node2, node3, node0);
	if (angle < min_angle) min_angle = angle;
	if (angle > max_angle) max_angle = angle;

	return 1.0 - std::max((max_angle - M_PI_HALF)/(M_PI - M_PI_HALF), (M_PI_HALF - min_angle)/(M_PI_HALF));
}

double MeshQualityEquiAngleSkew::checkTetrahedron (CElem const * const elem) const
{
	double const * const node0 (elem->GetNode(0)->getData());
	double const * const node1 (elem->GetNode(1)->getData());
	double const * const node2 (elem->GetNode(2)->getData());
	double const * const node3 (elem->GetNode(3)->getData());

	// first triangle (0,1,2)
	double min_angle (MATHLIB::getAngle (node1, node0, node2));
	double max_angle (min_angle);

	double angle (MATHLIB::getAngle (node2, node1, node0));
	if (angle < min_angle) min_angle = angle;
	if (angle > max_angle) max_angle = angle;

	angle = MATHLIB::getAngle (node0, node2, node1);
	if (angle < min_angle) min_angle = angle;
	if (angle > max_angle) max_angle = angle;

	// second triangle (0,1,3)
	angle = MATHLIB::getAngle (node3, node0, node1);
	if (angle < min_angle) min_angle = angle;
	if (angle > max_angle) max_angle = angle;

	angle = MATHLIB::getAngle (node0, node1, node3);
	if (angle < min_angle) min_angle = angle;
	if (angle > max_angle) max_angle = angle;

	angle = MATHLIB::getAngle (node1, node3, node0);
	if (angle < min_angle) min_angle = angle;
	if (angle > max_angle) max_angle = angle;

	// third triangle (0,2,3)
	angle = MATHLIB::getAngle (node3, node0, node2);
	if (angle < min_angle) min_angle = angle;
	if (angle > max_angle) max_angle = angle;

	angle = MATHLIB::getAngle (node0, node2, node3);
	if (angle < min_angle) min_angle = angle;
	if (angle > max_angle) max_angle = angle;

	angle = MATHLIB::getAngle (node2, node3, node0);
	if (angle < min_angle) min_angle = angle;
	if (angle > max_angle) max_angle = angle;

	// fourth triangle (1,2,3)
	angle = MATHLIB::getAngle (node3, node1, node2);
	if (angle < min_angle) min_angle = angle;
	if (angle > max_angle) max_angle = angle;

	angle = MATHLIB::getAngle (node1, node2, node3);
	if (angle < min_angle) min_angle = angle;
	if (angle > max_angle) max_angle = angle;

	angle = MATHLIB::getAngle (node2, node3, node1);
	if (angle < min_angle) min_angle = angle;
	if (angle > max_angle) max_angle = angle;

	return 1.0 - std::max((max_angle - M_PI_HALF)/(M_PI - M_PI_THIRD), (M_PI_THIRD - min_angle)/(M_PI_THIRD));
}

double MeshQualityEquiAngleSkew::checkHexahedron (CElem const * const elem) const
{
	return 0.0;
}

} // end namespace Mesh_Group
