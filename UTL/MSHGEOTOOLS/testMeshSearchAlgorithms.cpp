/*
 * testMeshSearchAlgorithms.cpp
 *
 *  Created on: Feb 28, 2012
 *      Author: TF
 */

#include <ctime>

// FileIO
#include "MeshIO/OGSMeshIO.h"
#include "XmlIO/XmlGmlInterface.h"

// GeoLib
#include "OctTree.h"
#include "AxisAlignedBoundingBox.h"
#include "GEOObjects.h"
#include "ProjectData.h"

// MathLib
#include "AnalyticalGeometry.h"
#include "Vector3.h"

// MSH
#include "msh_lib.h" // for FEMRead

// we need this for using the xml functions of Qt
#include <QApplication>

void testOctTree(MeshLib::CFEMesh const*const mesh)
{
// get the mesh node vector
	std::vector<MeshLib::CNode*> const& nodes_oct_tree(mesh->getNodeVector());
//	std::vector<GEOLIB::Point*> nodes_as_pnts;
	const size_t n_nodes_oct_tree(nodes_oct_tree.size());
	MeshLib::CNode min(0, nodes_oct_tree[0]->getData()), max(1, nodes_oct_tree[0]->getData());
//	GEOLIB::Point min(nodes_oct_tree[0]->getData()), max(nodes_oct_tree[0]->getData());
	// determine bounding box
	for (size_t k(0); k<n_nodes_oct_tree; k++) {
		if ((*nodes_oct_tree[k])[0] < min[0]) min[0] = (*nodes_oct_tree[k])[0];
		if ((*nodes_oct_tree[k])[1] < min[1]) min[1] = (*nodes_oct_tree[k])[1];
		if ((*nodes_oct_tree[k])[2] < min[2]) min[2] = (*nodes_oct_tree[k])[2];

		if (max[0] < (*nodes_oct_tree[k])[0]) max[0] = (*nodes_oct_tree[k])[0];
		if (max[1] < (*nodes_oct_tree[k])[1]) max[1] = (*nodes_oct_tree[k])[1];
		if (max[2] < (*nodes_oct_tree[k])[2]) max[2] = (*nodes_oct_tree[k])[2];

//		nodes_as_pnts.push_back(new GEOLIB::Point(nodes_oct_tree[k]->getData()));
	}
	GEOLIB::OctTree<MeshLib::CNode> oct_tree(min, max, 2);
	std::cout << "[OctTree] inserting " << n_nodes_oct_tree << " points ... " << std::flush;
	clock_t start(clock());
//	GEOLIB::OctTree<GEOLIB::Point> oct_tree(min, max, 2);
	for (size_t k(0); k<n_nodes_oct_tree; k++) {
		oct_tree.addPoint(nodes_oct_tree[k]);
//		oct_tree.addPoint(nodes_as_pnts[k]);
	}
	clock_t stop(clock());
	std::cout << "done,  " << (stop-start)/(double)(CLOCKS_PER_SEC) << " seconds" << std::endl;

	std::vector<MeshLib::CNode*> pnts_in_range;
	MeshLib::CNode q_min(min);
	MeshLib::CNode q_max(max);
//	std::vector<GEOLIB::Point*> pnts_in_range;
//	GEOLIB::Point q_min(min);
//	GEOLIB::Point q_max(max);
	start = clock();
	oct_tree.getPointsInRange(q_min, q_max, pnts_in_range);
	stop = clock();
	std::cout << "1st query: " << pnts_in_range.size() << " points in range " << q_min << " x " << q_max << " took " << (stop-start)/(double)(CLOCKS_PER_SEC) << " seconds" <<std::endl;

	pnts_in_range.clear();

	q_min[0] = (max[0] + min[0]) / 2.0;
	q_min[1] = (max[1] + min[1]) / 2.0;
	q_min[2] = (max[2] + min[2]) / 2.0;
	oct_tree.getPointsInRange(q_min, q_max, pnts_in_range);
	stop = clock();
	std::cout << "2nd query: " << pnts_in_range.size() << " points in range " << q_min << " x " << q_max << " took " << (stop-start)/(double)(CLOCKS_PER_SEC) << " seconds" <<std::endl;

	q_min[0] = min[0];
	q_min[1] = min[1];
	q_min[2] = min[2];
	q_max[0] = (max[0] + min[0]) / 2.0;
	q_max[1] = (max[1] + min[1]) / 2.0;
	q_max[2] = (max[2] + min[2]) / 2.0;
	start = clock();
	oct_tree.getPointsInRange(q_min, q_max, pnts_in_range);
	stop = clock();
	std::cout << "3rd query: " << pnts_in_range.size() << " points in range " << q_min << " x " << q_max << " took " << (stop-start)/(double)(CLOCKS_PER_SEC) << " seconds" <<std::endl;
}

int main (int argc, char* argv[])
{
	// Creating a non-gui (console) Qt application
	QApplication app(argc, argv, false);

	if (argc < 6) {
		std::cout << "Usage: " << argv[0] << " --mesh-in meshfile "
						<< std::flush;
		std::cout << "--geo-in geo-file --geo-out geo-out-file"
						<< std::endl;
		return -1;
	}

	// *** parsing mesh options from command line
	// *** parsing mesh in file name
	std::string tmp(argv[1]);
	if (tmp.find("--mesh-in") == std::string::npos) {
		std::cout << "could not find switch for reading mesh file name" << std::endl;
		return -1;
	}
	tmp = argv[2];
	std::string file_base_name(tmp);
	if (tmp.find(".msh") != std::string::npos) file_base_name = tmp.substr(0, tmp.size() - 4);

	// *** parsing geometry options from command line
	tmp = argv[3];
	if (tmp.find("--geo-in") == std::string::npos) {
		std::cout << "could not find switch for reading geometry file name" << std::endl;
		return -1;
	}
	std::string geo_fname_in(argv[4]);

	// *** parsing geo output name
	tmp = argv[5];
	if (tmp.find("--geo-out") == std::string::npos) {
		std::cout << "could not find switch for file name for writing the geometry" << std::endl;
		return -1;
	}
	std::string geo_fname_out(argv[6]);

	// *** read mesh
	std::vector<MeshLib::CFEMesh*> mesh_vec;
	FEMRead(file_base_name, mesh_vec);
	if (mesh_vec.empty()) {
		std::cerr << "could not read mesh from file " << std::endl;
		return -1;
	}
	MeshLib::CFEMesh* mesh (mesh_vec[mesh_vec.size() - 1]);
	mesh->setNumberOfNodesFromNodesVectorSize();

	// *** read geometry
	GEOLIB::GEOObjects* geo_objs (new GEOLIB::GEOObjects);
	ProjectData* project_data (new ProjectData);
	project_data->setGEOObjects (geo_objs);

	std::string schema_name("./OpenGeoSysGLI.xsd");
	FileIO::XmlGmlInterface xml(project_data, schema_name);
	xml.readFile(QString::fromStdString (geo_fname_in));
	std::vector<std::string> original_geo_names;
	geo_objs->getGeometryNames(original_geo_names);

	// *** begin test OctTree
	testOctTree(mesh);

//	xml.writeFile(geo_fname_out, original_geo_names[1]);

	delete project_data;

	return 0;
}
