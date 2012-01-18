/*
 * OGSMeshIO.h
 *
 *  Created on: Mar 3, 2011
 *      Author: TF
 */

#ifndef OGSMESHIO_H_
#define OGSMESHIO_H_

#include <sstream>
#include <iostream>

namespace MeshLib
{
class CFEMesh;
class CElem;
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
	
	/// Write an OGS msh to a stringstream
	static void write (MeshLib::CFEMesh const* mesh, std::stringstream &out);
private:
	static void writeElementsExceptLines (std::vector<MeshLib::CElem*> const& ele_vec, std::stringstream &out);
};
} // end namespace FileIO

#endif /* OGSMESHIO_H_ */
