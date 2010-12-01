/**
 * \file ProjectData.h
 * 25/08/2010 KR Initial implementation
 */

#ifndef PROJECTDATA_H_
#define PROJECTDATA_H_

#include "GEOObjects.h"
#include "msh_mesh.h"



class ProjectData
{
public:
	ProjectData();
	virtual ~ProjectData();

	// Returns the GEOObjects containing all points, polylines and surfaces
	GEOLIB::GEOObjects& getGEOObjects() const;

	/// Adds a new mesh
	virtual void addMesh(Mesh_Group::CFEMesh* mesh, std::string &name);

	/// Returns the mesh with the given name.
	const Mesh_Group::CFEMesh* getMesh(const std::string &name) const;

	/// Removes the mesh with the given name.
	virtual bool removeMesh(const std::string &name);

	/// Checks if the name of the mesh is already exists, if so it generates a unique name.
	bool isUniqueMeshName(std::string &name);

private:
	GEOLIB::GEOObjects _geoObjects;
	std::map<std::string, Mesh_Group::CFEMesh*> _msh_vec;
};

#endif //PROJECTDATA_H_
