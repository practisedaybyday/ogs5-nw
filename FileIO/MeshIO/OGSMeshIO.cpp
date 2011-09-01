/*
 * OGSMeshIO.cpp
 *
 *  Created on: Mar 3, 2011
 */

#include "MeshIO/OGSMeshIO.h"
#include "GEOObjects.h"
#include "msh_mesh.h"
#include "msh_lib.h"

namespace FileIO {

MeshLib::CFEMesh* OGSMeshIO::loadMeshFromFile(std::string const& fname)
{
	std::cout << "FEMRead ... " << std::flush;
/*
#ifndef NDEBUG
	QTime myTimer;
	myTimer.start();
#endif
*/
	FEMDeleteAll();

	MeshLib::CFEMesh* msh = FEMRead(fname.substr(0, fname.length()-4));
	if (msh)
	{
/*
#ifndef NDEBUG
		QTime constructTimer;
		constructTimer.start();
#endif
*/
		msh->ConstructGrid();

		std::cout << "Nr. Nodes: " << msh->nod_vector.size() << std::endl;
		std::cout << "Nr. Elements: " << msh->ele_vector.size() << std::endl;
/*
#ifndef NDEBUG
		std::cout << "constructGrid time: " << constructTimer.elapsed() << " ms" << std::endl;
#endif
*/
		msh->FillTransformMatrix();
/*
#ifndef NDEBUG
		std::cout << "Loading time: " << myTimer.elapsed() << " ms" << std::endl;
#endif
*/
		return msh;
	}

    std::cout << "Failed to load the mesh file: " << fname << std::endl;
	return NULL;
}

void OGSMeshIO::write(MeshLib::CFEMesh const * mesh, std::ofstream &out)
{
	out << "#FEM_MSH" << std::endl;

	out << "$PCS_TYPE" << std::endl << "  " << mesh->pcs_name << std::endl;

	out << "$NODES" << std::endl << "  ";
	out << mesh->GetNodesNumber(false) << std::endl;
	for (size_t i(0); i < mesh->nod_vector.size(); i++) {
		out << i
			<< " " << mesh->nod_vector[i]->X()
			<< " " << mesh->nod_vector[i]->Y()
			<< " " << mesh->nod_vector[i]->Z() << std::endl;
	}

	out << "$ELEMENTS" << std::endl << "  ";
	const size_t ele_vector_size (mesh->ele_vector.size());
	const double epsilon (std::numeric_limits<double>::epsilon());
	std::vector<bool> non_null_element (ele_vector_size, true);
	size_t non_null_elements (0);
	for (size_t i(0); i < ele_vector_size; i++) {
		if ((mesh->ele_vector[i])->calcVolume() < epsilon) { // || (mesh->ele_vector[i])->GetElementType() == MshElemType::LINE) {
			non_null_element[i] = false;
		} else {
			non_null_elements++;
		}
	}
	out << non_null_elements << std::endl;
	for (size_t i(0), k(0); i < ele_vector_size; i++) {
		if (non_null_element[i]) {
			mesh->ele_vector[i]->SetIndex(k);
			mesh->ele_vector[i]->WriteIndex(out);
			k++;
		}
	}

	out << " $LAYER" << std::endl;
	out << "  ";
	out << mesh->_n_msh_layer << std::endl;
	out << "#STOP" << std::endl;
}

} // end namespace FileIO
