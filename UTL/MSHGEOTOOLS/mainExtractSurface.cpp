/*
 * mainExtractSurface.cpp
 *
 *  Created on: Jan 26, 2011
 *      Author: TF
 */

#include "ExtractSurface.h"

// FileIO
#include "OGSIOVer4.h"
#include "XMLInterface.h"

// GEO
#include "GEOObjects.h"
#include "PolylineVec.h"
#include "Polygon.h"

// MSH
#include "msh_mesh.h"
#include "msh_lib.h" // for FEMRead

// STL
#include <string>
Problem *aproblem = NULL;


int main (int argc, char *argv[])
{
	if (argc < 5) {
		std::cout << "Usage: " << argv[0] << " --mesh ogs_meshfile --geometry ogs_geometry" << std::endl;
		return -1;
	}

	// *** read mesh
	std::string tmp (argv[1]);
	if (tmp.find ("--mesh") == std::string::npos) {
		std::cout << "could not extract mesh file name" << std::endl;
		return -1;
	}

	tmp = argv[2];
	std::string file_base_name (tmp);
	if (tmp.find (".msh") != std::string::npos)
		file_base_name = tmp.substr (0, tmp.size()-4);

	Mesh_Group::CFEMesh* mesh (FEMRead(file_base_name));
	if (!mesh) {
		std::cerr << "could not read mesh from file " << std::endl;
		return -1;
	}

	// *** read geometry
	tmp = argv[3];
	if (tmp.find ("--geometry") == std::string::npos) {
		std::cout << "could not extract geometry file name" << std::endl;
		return -1;
	}

	GEOLIB::GEOObjects *geo (new GEOLIB::GEOObjects);
	tmp = argv[4];
	FileIO::readGLIFileV4(tmp, geo);

	// *** get Polygon
	const std::vector<GEOLIB::Polyline*>* plys (geo->getPolylineVec (tmp));
	if (!plys) {
		std::cout << "could not get vector of polylines" << std::endl;
		delete mesh;
		delete geo;
		return -1;
	}

	std::vector<GEOLIB::Point*> *sfc_pnts (new std::vector<GEOLIB::Point*>);
	std::vector<GEOLIB::Surface*> *surfaces (new std::vector<GEOLIB::Surface*>);

	Mesh_Group::ExtractSurface extract_surface (mesh);
	GEOLIB::Polygon polygon (*((*plys)[0]));
	GEOLIB::Surface* surface (extract_surface.extractSurface(polygon, *sfc_pnts));
	std::cout << "surface: " << surface << std::endl;
	surfaces->push_back (surface);
	std::cout << "number of triangles in surface: " <<  surface->getNTriangles() << std::endl;

	std::string fname ("PointsForSurface");
	geo->addPointVec (sfc_pnts, fname);
	geo->addSurfaceVec (surfaces, fname);
//	// remove initial geometry
//	geo->removeSurfaceVec (tmp);
//	geo->removePolylineVec (tmp);
//	geo->removePointVec (tmp);

	XMLInterface xml_out (geo, "OpenGeoSysGLI.xsd");
	std::string out_fname ("Surface.gml");
	xml_out.writeGLIFile (out_fname , fname);

	delete mesh;
	delete geo;
}
