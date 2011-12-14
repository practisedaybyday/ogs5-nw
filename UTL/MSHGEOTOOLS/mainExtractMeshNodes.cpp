/*
 * extractMeshNodes.cpp
 *
 *  Created on: Dec 6, 2010
 *      Author: TF
 *
 */

// Base
#include "StringTools.h"

// FEM
#include "problem.h"

// MSH
#include "msh_lib.h" // for FEMRead
#include "msh_mesh.h"

// MSHGEOTOOLS
#include "ExtractMeshNodes.h"

// GEO
#include "GEOObjects.h"
#include "Polygon.h"
#include "PolylineVec.h"

// FileIO
#include "OGSIOVer4.h"

Problem* aproblem = NULL;

// MathLib
#include "Matrix.h"
#include "Vector3.h"

void getMeshNodesFromLayerAlongPolyline(MeshLib::CFEMesh *const mesh, GEOLIB::GEOObjects* geo,
				std::string const& name, size_t ply_id, size_t layer, std::vector<size_t> &mesh_ids)
{
	// *** get polyline from geo
	std::vector<GEOLIB::Polyline*> const* polylines(geo->getPolylineVecObj(name)->getVector());
	if (polylines->size() <= ply_id) {
		std::cout << "the given polyline id is greater than the number of polylines ("
						<< polylines->size() << ")" << std::endl;
		return;
	}

	GEOLIB::Polyline const*const ply_in((*polylines)[ply_id]);
	GEOLIB::Polyline* ply_out(NULL);

	MeshLib::ExtractMeshNodes extract_mesh_nodes (mesh);
	extract_mesh_nodes.getProjectedPolylineFromPolyline(*ply_in, geo, name, ply_out, layer);
	mesh->GetNODOnPLY(ply_out, mesh_ids);

	delete ply_out;
}


void writeMeshNodes(MeshLib::CFEMesh const* const mesh, std::vector<size_t> const& node_ids,
				std::string const& fname_ids, std::string const& fname_gli, bool write_gli)
{
	std::ofstream os (fname_ids.c_str());
	if (!os) {
		std::cout << "could not open file " << fname_ids << " for writing mesh node ids" << std::endl;
		return;
	}

	const size_t node_ids_size (node_ids.size());
	for (size_t k(0); k<node_ids_size; k++) {
		os << node_ids[k] << std::endl;
	}
	os.close();

	if (write_gli) {
		os.open(fname_gli.c_str());
		if (!os) {
			std::cout << "could not open file " << fname_gli << " for writing gli points" << std::endl;
			return;
		}
		std::vector<MeshLib::CNode*> const& nodes (mesh->getNodeVector());
		os << "#POINTS" << std::endl;
		for (size_t k(0); k<node_ids_size; k++) {
			double const*const node(nodes[node_ids[k]]->getData());
			os << k << " " << node[0] << " " << node[1] << " " << node[2] << std::endl;
		}
		os << "#STOP" << std::endl;
		os.close();
	}
}

int main (int argc, char* argv[])
{
	if (argc < 5)
	{
		std::cout << "Usage: " << argv[0] <<
		" --mesh ogs_meshfile --geometry ogs_geometry" << std::endl;
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
	mesh->ConstructGrid();

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
	FileIO::readGLIFileV4(tmp, geo, unique_name);

//	{
//		const std::vector<GEOLIB::Point*>* pnts (geo->getPointVec (tmp));
//		if (pnts) {
//			std::string fname ("MeshIDs.txt");
//			std::ofstream out (fname.c_str());
//
//			std::string fname_gli ("MeshNodesAsPnts.gli");
//			std::ofstream pnt_out (fname_gli.c_str());
//			pnt_out << "#POINTS" << std::endl;
//
//			MeshLib::ExtractMeshNodes extract_mesh_nodes (mesh);
//
//			const size_t n_pnts (pnts->size());
//			for (size_t k(0); k<n_pnts; k++) {
//				extract_mesh_nodes.writeNearestMeshNodeToPoint (out, pnt_out, *((*pnts)[k]));
//			}
//			pnt_out << "#STOP" << std::endl;
//		}
//		return 0;
//	}

	// *** get Polygon
	const std::vector<GEOLIB::Polyline*>* plys (geo->getPolylineVec (unique_name));
	if (!plys)
	{
		std::cout << "could not get vector of polylines" << std::endl;
		delete mesh;
		delete geo;
		return -1;
	}

	std::vector<size_t> mesh_ids;
	size_t ply_id (0);
	size_t layer(1);
	getMeshNodesFromLayerAlongPolyline(mesh, geo, unique_name, ply_id, layer, mesh_ids);
	writeMeshNodes(mesh, mesh_ids, "MeshIDs.txt", "MeshNodesAsPoints.gli", true);

	//*** extract surface out of mesh
//	MeshLib::ExtractMeshNodes extract_mesh_nodes (mesh);
//
//	// *** generate a surface from polyline
//	std::vector<GEOLIB::Polyline*> polylines;
//	const size_t n_plys (plys->size());
//	for (size_t k(0); k < n_plys; k++)
//	{
//		bool closed ((*plys)[k]->isClosed());
//		if (!closed && k == 19)
//		{
//			std::cout << "converting polyline " << k << " to closed polyline" <<
//			std::endl;
//			GEOLIB::Polygon* polygon (NULL);
//			extract_mesh_nodes.getPolygonFromPolyline(*((*plys)[k]), geo, tmp, polygon);
//			polylines.push_back (polygon);
//		}
//	}
//
//	geo->appendPolylineVec (polylines, tmp);
//	FileIO::writeGLIFileV4 ("New.gli", tmp, *geo);

	// *** search mesh nodes for direct assigning bc, st or ic
//	std::string fname ("MeshIDs.txt");
//	std::ofstream out (fname.c_str());
//
//	std::string fname_gli ("MeshNodesAsPnts.gli");
//	std::ofstream pnt_out (fname_gli.c_str());
//	pnt_out << "#POINTS" << std::endl;
//
//	const size_t n_plys (plys->size());
//	for (size_t k(0); k<n_plys; k++) {
//		bool closed ((*plys)[k]->isClosed());
//		if (!closed) {
//			std::cout << "polyline " << k << " is not closed" << std::endl;
//		} else {
//			GEOLIB::Polygon polygon (*((*plys)[k]));
////			extract_mesh_nodes.writeMesh2DNodeIDAndArea (out, pnt_out, polygon);
//			extract_mesh_nodes.writeTopSurfaceMeshNodeIDs (out, pnt_out, polygon);
//			// write all nodes - not only the surface nodes
////			extract_mesh_nodes.writeMeshNodeIDs (out, pnt_out, polygon);
//		}
//	}

	// *** for Model Pipiripau
//	std::vector<GEOLIB::Polygon*> holes;
//	size_t bounding_polygon_id(0);
//	while (bounding_polygon_id < n_plys && ! (*plys)[bounding_polygon_id]->isClosed()) {
//		bounding_polygon_id++;
//	}
//
//	for (size_t k(bounding_polygon_id+1); k<n_plys; k++) {
//		bool closed ((*plys)[k]->isClosed());
//		if (!closed) {
//			std::cout << "polyline " << k << " is not closed" << std::endl;
//		} else {
//			holes.push_back (new GEOLIB::Polygon(*(((*plys)[k]))));
//		}
//	}
//	extract_mesh_nodes.writeMesh2DNodeIDAndArea (out, pnt_out, GEOLIB::Polygon((*((*plys)[bounding_polygon_id]))), holes);
//	for (size_t k(0); k<holes.size(); k++) {
//		delete holes[k];
//	}

//	out << "#STOP" << std::endl;
//	out.close();
//
//	pnt_out << "#STOP" << std::endl;
//	pnt_out.close ();

	delete mesh;
	delete geo;

	return 0;
}
