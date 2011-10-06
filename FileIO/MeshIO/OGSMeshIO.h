/*
 * OGSMeshIO.h
 *
 *  Created on: Mar 3, 2011
 *      Author: TF
 */

#ifndef OGSMESHIO_H_
#define OGSMESHIO_H_

#include <fstream>

namespace MeshLib
{
class CFEMesh;
}

namespace FileIO
{
class OGSMeshIO
{
public:
	/// Read a OGS mesh from file.
	static MeshLib::CFEMesh* loadMeshFromFile(std::string const& fileName);

	/// Write an OGS msh to a file.
	static void write (MeshLib::CFEMesh const* mesh, std::ofstream &out);
};
} // end namespace FileIO

#endif /* OGSMESHIO_H_ */
