/*
 * MeshQualityChecker.h
 *
 *  Created on: Dec 8, 2010
 *      Author: TF
 */

#ifndef MESHQUALITYCHECKER_H_
#define MESHQUALITYCHECKER_H_

#include <vector>

// MSH
#include "msh_mesh.h"

namespace Mesh_Group {

class MeshQualityChecker {
public:
	MeshQualityChecker(CFEMesh const * const mesh);
	~MeshQualityChecker();

	void checkTriangles ();

	const std::vector<size_t>& getHistogramm () const { return _static_histogramm; }

private:
	CFEMesh const * const _mesh;
	std::vector<size_t> _static_histogramm;
};

}

#endif /* MESHQUALITYCHECKER_H_ */
