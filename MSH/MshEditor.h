/**
 * \file MshEditor.h
 * 2011/06/15 KR Initial implementation
 */

#ifndef MSHEDITOR_H
#define MSHEDITOR_H

#include <cstddef>
#include <vector>

namespace GEOLIB {
	class PointWithID;
}

namespace MeshLib
{
class CFEMesh;
}

class GridAdapter;

/**
 * \brief A set of tools for manipulating existing meshes
 */
class MshEditor
{
public:
	MshEditor() {}
	~MshEditor() {}

	static MeshLib::CFEMesh* removeMeshNodes(MeshLib::CFEMesh* mesh,
	                                         const std::vector<size_t> &nodes);

	static const std::vector<GEOLIB::PointWithID*> getSurfaceNodes(const MeshLib::CFEMesh &mesh);

	static MeshLib::CFEMesh* getMeshSurface(const MeshLib::CFEMesh &mesh);

private:
};

#endif //MSHEDITOR_H
