/*
 * OGSmsh2GeoIO.h
 *
 *  Created on: Sep 3, 2012
 *      Author: NW
 */

#ifndef OGSMESH2GEOIO_H_
#define OGSMESH2GEOIO_H_

#include <string>

namespace GEOLIB
{
class GEOObjects;
}

namespace MeshLib
{
class CFEMesh;
}

namespace FileIO
{
class OGSmsh2GeoIO
{
public:
	/**
	 * load a surface mesh (ogs mesh format) as a geometric surface
     * @param mesh OGS mesh object
	 * @param fname file name of the surface mesh
	 * @param geo the object that manages all geometries,
	 * new surface will be put into this container
	 */
	static void loadMeshAsGeometry (const MeshLib::CFEMesh* mesh, std::string &project_name,GEOLIB::GEOObjects* geo);
};
} // end namespace FileIO

#endif /* OGSMESH2GEOIO_H_ */
