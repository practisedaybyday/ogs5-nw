/*
 * OGSMeshIO.h
 *
 *  Created on: Mar 3, 2011
 *      Author: TF
 */

#ifndef OGSMESHIO_H_
#define OGSMESHIO_H_

#include "Writer.h"

#include <sstream>
#include <iostream>
#include <vector>

namespace MeshLib
{
class CFEMesh;
class CElem;
}

namespace FileIO
{
class OGSMeshIO : public Writer
{
public:
	/// @brief Constructor.
	OGSMeshIO();

	/// @brief Read a OGS mesh from file.
	MeshLib::CFEMesh* loadMeshFromFile(std::string const& fileName);

	/// @brief Write functionality.
	void write(std::ostream &out);

	/// @brief Sets the mesh.
	void setMesh(MeshLib::CFEMesh const* mesh);

private:
	void writeElementsExceptLines (std::vector<MeshLib::CElem*> const& ele_vec, std::ostream &out);

	MeshLib::CFEMesh const* _mesh;
};
} // end namespace FileIO

#endif /* OGSMESHIO_H_ */
