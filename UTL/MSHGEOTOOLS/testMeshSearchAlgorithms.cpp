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

void testOctTreeGlobalQueries (GEOLIB::OctTree<MeshLib::CNode> const& oct_tree, double const*const min, double const*const max)
{
	std::vector<MeshLib::CNode*> pnts_in_range;
	MeshLib::CNode q_min(min);
	MeshLib::CNode q_max(max);
//	std::vector<GEOLIB::Point*> pnts_in_range;
//	GEOLIB::Point q_min(min);
//	GEOLIB::Point q_max(max);
	clock_t start, stop;

	start = clock();
	oct_tree.getPointsInRange(q_min, q_max, pnts_in_range);
	stop = clock();
	std::cout << "query all: " << pnts_in_range.size() << " points in range " << q_min << " x " << q_max << " took " << (stop-start)/(double)(CLOCKS_PER_SEC) << " seconds" <<std::endl;
	pnts_in_range.clear();

	size_t cnt(0);
	q_min[0] = (max[0] + min[0]) / 2.0;
	q_min[1] = (max[1] + min[1]) / 2.0;
	q_min[2] = (max[2] + min[2]) / 2.0;
	start = clock();
	oct_tree.getPointsInRange(q_min, q_max, pnts_in_range);
	stop = clock();
	std::cout << "query NEU: " << pnts_in_range.size() << " points in range " << q_min << " x " << q_max << " took " << (stop-start)/(double)(CLOCKS_PER_SEC) << " seconds" <<std::endl;
	cnt += pnts_in_range.size();
	pnts_in_range.clear();

	q_min[0] = min[0];
	q_min[1] = (max[1] + min[1]) / 2.0;
	q_min[2] = (max[2] + min[2]) / 2.0;
	q_max[0] = (max[0] + min[0]) / 2.0;
	q_max[1] = max[1];
	q_max[2] = max[2];
	start = clock();
	oct_tree.getPointsInRange(q_min, q_max, pnts_in_range);
	stop = clock();
	std::cout << "query NWU: " << pnts_in_range.size() << " points in range " << q_min << " x " << q_max << " took " << (stop-start)/(double)(CLOCKS_PER_SEC) << " seconds" <<std::endl;
	cnt += pnts_in_range.size();
	pnts_in_range.clear();

	q_min[0] = min[0];
	q_min[1] = min[1];
	q_min[2] = (min[2]+max[2])/2.0;
	q_max[0] = (max[0] + min[0]) / 2.0;
	q_max[1] = (max[1] + min[1]) / 2.0;
	q_max[2] = max[2];
	start = clock();
	oct_tree.getPointsInRange(q_min, q_max, pnts_in_range);
	stop = clock();
	std::cout << "query SWU: " << pnts_in_range.size() << " points in range " << q_min << " x " << q_max << " took " << (stop-start)/(double)(CLOCKS_PER_SEC) << " seconds" <<std::endl;
	cnt += pnts_in_range.size();
	pnts_in_range.clear();

	q_min[0] = (max[0]+min[0])/2.0;
	q_min[1] = min[1];
	q_min[2] = (min[2]+max[2])/2.0;
	q_max[0] = max[0];
	q_max[1] = (max[1] + min[1]) / 2.0;
	q_max[2] = max[2];
	start = clock();
	oct_tree.getPointsInRange(q_min, q_max, pnts_in_range);
	stop = clock();
	std::cout << "query SEU: " << pnts_in_range.size() << " points in range " << q_min << " x " << q_max << " took " << (stop-start)/(double)(CLOCKS_PER_SEC) << " seconds" <<std::endl;
	cnt += pnts_in_range.size();
	pnts_in_range.clear();


	q_min[0] = (max[0] + min[0]) / 2.0;
	q_min[1] = (max[1] + min[1]) / 2.0;
	q_min[2] = min[2];
	q_max[0] = max[0];
	q_max[1] = max[1];
	q_max[2] = (max[2]+min[2])/2.0;
	start = clock();
	oct_tree.getPointsInRange(q_min, q_max, pnts_in_range);
	stop = clock();
	std::cout << "query NEL: " << pnts_in_range.size() << " points in range " << q_min << " x " << q_max << " took " << (stop-start)/(double)(CLOCKS_PER_SEC) << " seconds" <<std::endl;
	cnt += pnts_in_range.size();
	pnts_in_range.clear();

	q_min[0] = min[0];
	q_min[1] = (max[1] + min[1]) / 2.0;
	q_min[2] = min[2];
	q_max[0] = (max[0] + min[0]) / 2.0;
	q_max[1] = max[1];
	q_max[2] = (max[2]+min[2])/2.0;
	start = clock();
	oct_tree.getPointsInRange(q_min, q_max, pnts_in_range);
	stop = clock();
	std::cout << "query NWL: " << pnts_in_range.size() << " points in range " << q_min << " x " << q_max << " took " << (stop-start)/(double)(CLOCKS_PER_SEC) << " seconds" <<std::endl;
	cnt += pnts_in_range.size();
	pnts_in_range.clear();

	q_min[0] = min[0];
	q_min[1] = min[1];
	q_min[2] = min[2];
	q_max[0] = (max[0] + min[0]) / 2.0;
	q_max[1] = (max[1] + min[1]) / 2.0;
	q_max[2] = (max[2]+min[2])/2.0;
	start = clock();
	oct_tree.getPointsInRange(q_min, q_max, pnts_in_range);
	stop = clock();
	std::cout << "query SWL: " << pnts_in_range.size() << " points in range " << q_min << " x " << q_max << " took " << (stop-start)/(double)(CLOCKS_PER_SEC) << " seconds" <<std::endl;
	cnt += pnts_in_range.size();
	pnts_in_range.clear();

	q_min[0] = (max[0]+min[0])/2.0;
	q_min[1] = min[1];
	q_min[2] = min[2];
	q_max[0] = max[0];
	q_max[1] = (max[1] + min[1]) / 2.0;
	q_max[2] = (max[2]+min[2])/2.0;
	start = clock();
	oct_tree.getPointsInRange(q_min, q_max, pnts_in_range);
	stop = clock();
	std::cout << "query SEL: " << pnts_in_range.size() << " points in range " << q_min << " x " << q_max << " took " << (stop-start)/(double)(CLOCKS_PER_SEC) << " seconds" <<std::endl;
	cnt += pnts_in_range.size();
	pnts_in_range.clear();

	std::cout << cnt << " points in children" << std::endl;
}

void testOctTreePointQueries(GEOLIB::OctTree<MeshLib::CNode> const& oct_tree,
				std::vector<GEOLIB::Point*> const& pnts_for_search,
				std::vector<size_t> &idx_found_nodes)
{
	const size_t n_pnts(pnts_for_search.size());
	std::vector<MeshLib::CNode*> pnts_in_range;
	for(size_t k(0); k<n_pnts; k++) {
		MeshLib::CNode q_min (pnts_for_search[k]->getData());
		MeshLib::CNode q_max (pnts_for_search[k]->getData());
		for (size_t d(0); d<3; d++) {
			q_min[d] -= 1; //fabs(q_min[d]) * std::numeric_limits<double>::epsilon();
			q_max[d] += 1; // fabs(q_max[d]) * std::numeric_limits<double>::epsilon();
		}
		oct_tree.getPointsInRange(q_min, q_max, pnts_in_range);
		if (! pnts_in_range.empty()) {
			const size_t n_pnts_in_range(pnts_in_range.size());
			double sqr_dist(0.0);
			double distmin(MathLib::sqrDist (pnts_in_range[0]->getData(), pnts_for_search[k]->getData()));
			size_t idx(pnts_in_range[0]->GetIndex());
			for (size_t j(1); j<n_pnts_in_range; j++) {
				sqr_dist = MathLib::sqrDist (pnts_in_range[j]->getData(), pnts_for_search[k]->getData());
				if (sqr_dist < distmin) {
					distmin = sqr_dist;
					idx = pnts_in_range[j]->GetIndex();
				}
			}
			idx_found_nodes.push_back(idx);
		} else {
			std::cout << "no node found for point " << *(pnts_for_search[k]) << std::endl;
		}
		pnts_in_range.clear();
	}
}

void testOctTree(MeshLib::CFEMesh const*const mesh, std::vector<GEOLIB::Point*>& pnts_for_search, std::vector<size_t> &idx_found_nodes)
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

	// *** create OctTree object
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

//	testOctTreeGlobalQueries(oct_tree, min.getData(), max.getData());

	std::cout << "[OctTree] searching " << pnts_for_search.size() << " points ... " << std::flush;
	start = clock();
	testOctTreePointQueries(oct_tree, pnts_for_search, idx_found_nodes);
	stop = clock();
	std::cout << "done,  " << (stop-start)/(double)(CLOCKS_PER_SEC) << " seconds" << std::endl;
}

void testNaiveAlgorithm (MeshLib::CFEMesh const*const mesh,
				std::vector<GEOLIB::Point*>& pnts_for_search,
				std::vector<size_t> &idx_found_nodes)
{
	const size_t nodes_in_usage(static_cast<size_t> (mesh->NodesInUsage()));
	std::vector<MeshLib::CNode*> const& nodes (mesh->getNodeVector());
	const size_t n_pnts_for_search(pnts_for_search.size());
	std::cout << "[LinearSearchAlgorithm] searching " << pnts_for_search.size() << " points ... " << std::flush;
	clock_t start = clock();
	for (size_t k(0); k<n_pnts_for_search; k++) {
		double const*const k_th_pnt (pnts_for_search[k]->getData());
		double sqr_dist(0.0), distmin(MathLib::sqrDist (nodes[0]->getData(), k_th_pnt));
		size_t idx(0);
		for (size_t i = 1; i < nodes_in_usage; i++) {
			sqr_dist = MathLib::sqrDist (nodes[i]->getData(), k_th_pnt);
			if (sqr_dist < distmin) {
				distmin = sqr_dist;
				idx = i;
			}
		}
		idx_found_nodes.push_back(idx);
	}
	clock_t stop = clock();
	std::cout << "done,  " << (stop-start)/(double)(CLOCKS_PER_SEC) << " seconds" << std::endl;
}

void testMeshGridAlgorithm(MeshLib::CFEMesh const*const mesh,
				std::vector<GEOLIB::Point*>& pnts_for_search,
				std::vector<size_t> &idx_found_nodes)
{
	const size_t n_pnts_for_search(pnts_for_search.size());
	std::cout << "[MeshGridAlgorithm] searching " << pnts_for_search.size() << " points ... " << std::flush;
	clock_t start = clock();
	for (size_t k(0); k<n_pnts_for_search; k++) {
		idx_found_nodes.push_back(static_cast<size_t>(mesh->GetNODOnPNT(pnts_for_search[k])));
	}
	clock_t stop = clock();
	std::cout << "done,  " << (stop-start)/(double)(CLOCKS_PER_SEC) << " seconds" << std::endl;
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

	// *** preparing test data
	std::vector<MeshLib::CNode*> const& nodes(mesh->getNodeVector());
	std::vector<GEOLIB::Point*> pnts_for_search;
	const size_t n(10000);
	for (size_t k(0); k<n; k++) {
		double const*const c(nodes[2*k]->getData());
		// perturb the coordinates a little bit
		pnts_for_search.push_back(new GEOLIB::Point(c[0]+0.1, c[1]+0.5, c[2]+1.0));
	}

	// *** test mesh grid algorithm
	mesh->ConstructGrid();
	std::vector<size_t> idx_found_nodes_mesh_grid_alg;
	testMeshGridAlgorithm(mesh, pnts_for_search, idx_found_nodes_mesh_grid_alg);

	// *** test OctTree
	std::vector<size_t> idx_found_nodes_oct_tree;
	testOctTree(mesh, pnts_for_search, idx_found_nodes_oct_tree);

	// *** test linear algorithm
	std::vector<size_t> idx_found_nodes_linear_alg;
	testNaiveAlgorithm(mesh, pnts_for_search, idx_found_nodes_linear_alg);

	// *** compare results
	std::cout << "compare results ... " << std::flush;
	for (size_t k(0); k<idx_found_nodes_oct_tree.size(); k++) {
		if (idx_found_nodes_oct_tree[k] != idx_found_nodes_linear_alg[k]) {
			std::cout << std::endl << "point: " << *pnts_for_search[k] << " node found within oct_tree: "
				<<  *(nodes[idx_found_nodes_oct_tree[k]]) << ", node found with linear algorithm "
				<<  *(nodes[idx_found_nodes_linear_alg[k]]);
		}
	}
	std::cout << " done" << std::endl;

	for (size_t k(0); k<n; k++) {
		delete pnts_for_search[k];
	}

//	xml.writeFile(geo_fname_out, original_geo_names[1]);

	delete project_data;

	return 0;
}