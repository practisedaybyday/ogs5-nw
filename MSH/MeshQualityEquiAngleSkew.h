/*
 * MeshQualityEquiAngleSkew.h
 *
 *  Created on: Mar 17, 2011
 *      Author: TF
 */

#ifndef MESHQUALITYEQUIANGLESKEW_H_
#define MESHQUALITYEQUIANGLESKEW_H_

#include "MeshQualityChecker.h"

namespace Mesh_Group {

class MeshQualityEquiAngleSkew: public Mesh_Group::MeshQualityChecker {
public:
	MeshQualityEquiAngleSkew(CFEMesh const * const mesh);
	virtual ~MeshQualityEquiAngleSkew();

	virtual void check ();

private:
	double checkTriangle (CElem const * const elem) const;
	double checkQuad (CElem const * const elem) const;
	double checkTetrahedron (CElem const * const elem) const;
	double checkHexahedron (CElem const * const elem) const;

	const double M_PI_THIRD;
	const double M_PI_HALF;
};

}

#endif /* MESHQUALITYEQUIANGLESKEW_H_ */
