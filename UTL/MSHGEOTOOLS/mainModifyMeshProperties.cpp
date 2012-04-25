/*
 * modifyMeshProperties.cpp
 *
 *  Created on: Dec 8, 2010
 *      Author: TF
 *
 */

#include "ModifyMeshProperties.h"

// FEM
#include "problem.h"

// MSH
#include "msh_lib.h" // for FEMRead
#include "msh_mesh.h"

// GEO
#include "GEOObjects.h"
#include "Polygon.h"
#include "PolylineVec.h"

// FileIO
#include "MeshIO/OGSMeshIO.h"
#include "OGSIOVer4.h"

Problem* aproblem = NULL;

int main (int argc, char* argv[])
{
	if (argc < 9)
	{
		std::cout << "Usage: " << argv[0] <<
		" --mesh ogs_meshfile --geometry ogs_geometry --polygon_id id --material_id id" << std::endl;
		return -1;
	}

	// *** read mesh
	std::string tmp (argv[1]);
	if (tmp.find ("--mesh") == std::string::npos)
	{
		std::cout << "could not extract mesh file name" << std::endl;
		return -1;
	}

	tmp = argv[2];
	std::string file_base_name (tmp);
	if (tmp.find (".msh") != std::string::npos)
		file_base_name = tmp.substr (0, tmp.size() - 4);

	std::vector<MeshLib::CFEMesh*> mesh_vec;
	FEMRead(file_base_name, mesh_vec);
	if (mesh_vec.empty())
	{
		std::cerr << "could not read mesh from file " << std::endl;
		return -1;
	}
	MeshLib::CFEMesh* mesh (mesh_vec[mesh_vec.size() - 1]);

	// *** read geometry
	tmp = argv[3];
	if (tmp.find ("--geometry") == std::string::npos)
	{
		std::cout << "could not extract geometry file name" << std::endl;
		return -1;
	}

	GEOLIB::GEOObjects* geo (new GEOLIB::GEOObjects);
	tmp = argv[4];
	std::string unique_name;
	std::vector<std::string> error_strings;
	FileIO::readGLIFileV4(tmp, geo, unique_name, error_strings);

	// *** get Polygon
	const std::vector<GEOLIB::Polyline*>* plys (geo->getPolylineVec (tmp));
	if (!plys)
	{
		std::cout << "could not get vector of polylines" << std::endl;
		delete mesh;
		delete geo;
		return -1;
	}

	bool closed ((*plys)[0]->isClosed());
	if (!closed)
	{
		std::cout << "polyline is not closed" << std::endl;
		delete mesh;
		delete geo;
		return -1;
	}

	// *** get polygon id
	tmp = argv[5];
	if (tmp.find ("--polygon_id") == std::string::npos)
	{
		std::cout << "could not read polygon id" << std::endl;
		delete mesh;
		delete geo;
		return -1;
	}
	size_t polygon_id (atoi (argv[6]));

	// *** get material id
	tmp = argv[7];
	if (tmp.find ("--material") == std::string::npos)
	{
		std::cout << "could not read material id" << std::endl;
		delete mesh;
		delete geo;
		return -1;
	}
	size_t material_id (atoi (argv[8]));

	GEOLIB::Polygon polygon (*((*plys)[polygon_id]));
	MeshLib::ModifyMeshProperties modify_mesh_nodes (mesh);

	modify_mesh_nodes.setMaterial (polygon, material_id);

	std::string mesh_out_fname("MeshWithMaterial.msh");

	mesh->ConstructGrid();

	std::cout << "write " << mesh_out_fname << " ... " << std::flush;
	FileIO::OGSMeshIO mesh_io;
	mesh_io.setMesh(mesh);
	mesh_io.writeToFile (mesh_out_fname);
	std::cout << "ok" << std::endl;

	return 0;
}
