/*
 * MeshQualityNormalisedArea.h
 *
 * 2011/03/17 KR Initial Implementation
 */

#ifndef MESHQUALITYNORMALISEDAREA_H_
#define MESHQUALITYNORMALISEDAREA_H_

#include "MeshQualityChecker.h"

namespace Mesh_Group {

class MeshQualityNormalisedArea : public MeshQualityChecker {
public:
	MeshQualityNormalisedArea(CFEMesh const * const mesh);
	virtual ~MeshQualityNormalisedArea() {};

	virtual void check ();
};

}

#endif /* MESHQUALITYNORMALISEDAREA_H_ */
