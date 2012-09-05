/*
 * convertSurfaceMeshToGeometry.cpp
 *
 *  Created on: May 30, 2012
 *      Author: fischeth
 */

// STL
#include <vector>
#include <string>

// FileIO
#include "OGSIOVer4.h"
#include "XmlIO/XmlGmlInterface.h"

// FileIO/GeoIO
#include "GeoIO/OGSmsh2GeoIO.h"
// MSH
#include "msh_lib.h" // for FEMRead
#include "msh_mesh.h"

// GEO
#include "GEOObjects.h"
#include "Point.h"
#include "PolylineVec.h"
#include "ProjectData.h"

int main (int argc, char* argv[])
{
	if (argc < 5) {
		std::cout << "Usage: " << argv[0] << " --mesh ogs_mesh --geometry new_geo_file" << std::endl;
		return -1;
	}

	GEOLIB::GEOObjects* geo (new GEOLIB::GEOObjects);

	// *** read OGS mesh
    std::string tmp (argv[2]);
    std::string file_base_name (tmp);
    if (tmp.find (".msh") != std::string::npos)
        file_base_name = tmp.substr (0, tmp.size() - 4);
	std::vector<MeshLib::CFEMesh*> mesh_vec;
	FEMRead(file_base_name, mesh_vec);
	if (mesh_vec.empty())
	{
		std::cerr << "could not read mesh from file - " << file_base_name << std::endl;
		return -1;
	}
	MeshLib::CFEMesh* mesh (mesh_vec[mesh_vec.size() - 1]);
	if (mesh->getNodeVector().size() == 0 || mesh->getElementVector().size() == 0) {
        std::cerr << "could not read mesh from file - " << file_base_name << std::endl;
        return -1;
	}
	
	FileIO::OGSmsh2GeoIO::loadMeshAsGeometry(mesh, file_base_name, geo);

	ProjectData* project_data (new ProjectData);
	project_data->setGEOObjects (geo);
	FileIO::XmlGmlInterface xml_out (project_data, "OpenGeoSysGLI.xsd");

	std::string geo_fname(argv[4]);
	xml_out.setNameForExport(file_base_name);
	xml_out.writeToFile(geo_fname);

	// get Surface for writing a OpenGeoSys TIN
	std::vector<GEOLIB::Surface*> const& sfcs(*(geo->getSurfaceVec(file_base_name)));
	GEOLIB::Surface const*const sfc(sfcs[0]);
	std::ofstream out((geo_fname+".tin").c_str());
	for (size_t k(0); k<sfc->getNTriangles(); k++) {
		GEOLIB::Triangle const& tri_k(*((*sfc)[k]));
        out << k << " ";
		for (size_t j(0); j<3; j++) {
			out << *(tri_k.getPoint(j)) << " " << std::flush;
		}
		out << std::endl;
	}
	out.close();

	delete project_data;
}
