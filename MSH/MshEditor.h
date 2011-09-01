/**
 * \file MshEditor.h
 * 2011/06/15 KR Initial implementation
 */

#ifndef MSHEDITOR_H
#define MSHEDITOR_H

#include <vector>
#include <cstddef>

namespace MeshLib
{
	class CFEMesh;
}

/**
 * \brief A set of tools for manipulating existing meshes
 */
class MshEditor
{
public:
	MshEditor() {};
	~MshEditor() {};

	static MeshLib::CFEMesh* removeMeshNodes(MeshLib::CFEMesh* mesh, const std::vector<size_t> &nodes);

private:


};

#endif //MSHEDITOR_H
