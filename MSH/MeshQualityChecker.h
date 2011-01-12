/*
 * MeshQualityChecker.h
 *
 *  Created on: Dec 8, 2010
 *      Author: fischeth
 */

#ifndef MESHQUALITYCHECKER_H_
#define MESHQUALITYCHECKER_H_

#include <vector>

// MSH
#include "msh_mesh.h"

namespace Mesh_Group {

class MeshQualityChecker {
public:
	MeshQualityChecker(CFEMesh *mesh);
	~MeshQualityChecker();

	void checkTriangles ();

private:
	CFEMesh* _mesh;
	std::vector<size_t> _static_histogramm;
};

}

#endif /* MESHQUALITYCHECKER_H_ */
