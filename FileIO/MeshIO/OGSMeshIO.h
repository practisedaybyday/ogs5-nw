/*
 * OGSMeshIO.h
 *
 *  Created on: Mar 3, 2011
 *      Author: TF
 */

#ifndef OGSMESHIO_H_
#define OGSMESHIO_H_

#include <fstream>

namespace Mesh_Group {
class CFEMesh;
}

namespace FileIO {

class OGSMeshIO {
public:
	OGSMeshIO();
	virtual ~OGSMeshIO();

	void write (Mesh_Group::CFEMesh const * mesh, std::ofstream &out) const;
};

} // end namespace FileIO

#endif /* OGSMESHIO_H_ */
