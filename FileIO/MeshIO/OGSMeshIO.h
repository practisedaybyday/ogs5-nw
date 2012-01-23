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

	/// Write an OGS msh to a output stream.
	static void write (MeshLib::CFEMesh const* mesh, std::ostream &out);

private:
	static void writeElementsExceptLines (std::vector<MeshLib::CElem*> const& ele_vec, std::ostream &out);
};
} // end namespace FileIO

#endif /* OGSMESHIO_H_ */
