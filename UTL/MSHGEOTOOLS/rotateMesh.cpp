/*
 * rotateMesh.cpp
 *
 *  Created on: Jan 30, 2012
 *      Author: TF
 */

// FileIO
#include "MeshIO/OGSMeshIO.h"
#include "OGSIOVer4.h"

// MathLib
#include "AnalyticalGeometry.h"
#include "Vector3.h"

// MSH
#include "msh_lib.h" // for FEMRead

int main (int argc, char* argv[])
{
	if (argc < 5)
	{
		std::cout << "Usage: " << argv[0] <<
		" --mesh-in meshfile --mesh-out mesh-out-file" << std::endl;
		return -1;
	}

	// *** read mesh
	std::string tmp (argv[1]);
	if (tmp.find ("--mesh-in") == std::string::npos)
	{
		std::cout << "could not find switch for reading mesh file name" << std::endl;
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
	mesh->setNumberOfNodesFromNodesVectorSize();

	// *** read output name
	tmp = argv[3];
	if (tmp.find ("--mesh-out") == std::string::npos)
	{
		std::cout << "could not find switch for file name for writing the mesh" << std::endl;
		return -1;
	}

	// *** transfer mesh nodes to points
	// to make a copy is for this small programm useless, but it is more clear and
	// if this functions will become part of official OGS we are ready to use it
	MeshLib::CFEMesh mesh_copy(*mesh);

	std::vector<GEOLIB::Point*> pnts;
	const size_t n_nodes(mesh_copy.GetNodesNumber(false));
	std::vector<MeshLib::CNode*>& nodes(const_cast<std::vector<MeshLib::CNode*>& >(mesh_copy.getNodeVector()));
	for (size_t k(0); k<n_nodes; k++) {
		pnts.push_back(new GEOLIB::Point(nodes[k]->getData()));
	}

	MathLib::Vector plane_normal(0.0, 0.0, 0.0);
	double d(0.0);
	MathLib::getNewellPlane(pnts, plane_normal, d);
	MathLib::rotatePointsToXZ(plane_normal, pnts);

	double mean_val_y((*pnts[0])[1]);
	double min_val_y((*pnts[0])[1]), max_val_y((*pnts[0])[1]);
	for (size_t k(1); k<n_nodes; k++) {
		mean_val_y += (*pnts[k])[1];
		if ((*pnts[k])[1] < min_val_y) min_val_y = (*pnts[k])[1];
		if (max_val_y < (*pnts[k])[1]) max_val_y = (*pnts[k])[1];
	}
	mean_val_y /= n_nodes;

	double varianz (MathLib::fastpow(mean_val_y-(*pnts[0])[1],2));
	for (size_t k(1); k<n_nodes; k++) {
		varianz += MathLib::fastpow(mean_val_y-(*pnts[k])[1],2);
	}
	std::cout << "statistical data of y coordinate:" << std::endl;
	std::cout << "\tmean value: " << mean_val_y << std::endl;
	std::cout << "\tminimal value: " << min_val_y << std::endl;
	std::cout << "\tmaximal value: " << max_val_y << std::endl;
	std::cout << "\tvarianz: " << varianz << std::endl;
	std::cout << "\tstandard deviation: " << sqrt(varianz) << std::endl << std::endl;

	size_t answer(0);
	while (answer < 1 || answer > 3) {
		std::cout << "Would should I do? Please give the number of the alternative!" << std::endl;
		std::cout << "\t1 Leave the y coordinate for every mesh node untouched." << std::endl;
		std::cout << "\t2 Choose the mean value as y coordinate for every mesh node." << std::endl;
		std::cout << "\t3 Set y = 0 for every mesh node." << std::endl;
		std::cin >> answer;
	}

	if (answer == 1) {
		for (size_t k(0); k<n_nodes; k++) {
			nodes[k]->SetCoordinates(pnts[k]->getData());
		}
	} else {
		if (answer == 2) {
			for (size_t k(0); k<n_nodes; k++) {
				(*pnts[k])[1] = mean_val_y;
				nodes[k]->SetCoordinates(pnts[k]->getData());
			}
		} else {
			for (size_t k(0); k<n_nodes; k++) {
				(*pnts[k])[1] = 0.0;
				nodes[k]->SetCoordinates(pnts[k]->getData());
			}
		}
	}

	std::ofstream mesh_out;
	mesh_out.open (argv[4]);
	if (mesh_out.is_open())
	{
		std::cout << "writing rotated mesh to " << argv[4] << " ... " << std::flush;
		FileIO::OGSMeshIO::write (&mesh_copy, mesh_out);
		std::cout << "done" << std::endl;
	}
	mesh_out.close ();

	return 0;
}
