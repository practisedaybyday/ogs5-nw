/*
 * OGSMeshIO.cpp
 *
 *  Created on: Mar 3, 2011
 */

#include "MeshIO/OGSMeshIO.h"
#include "msh_mesh.h"

namespace FileIO {

OGSMeshIO::OGSMeshIO()
{}

OGSMeshIO::~OGSMeshIO()
{}

void OGSMeshIO::write(Mesh_Group::CFEMesh const * mesh, std::ofstream &out) const
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
	const double epsilon (sqrt(fabs(std::numeric_limits<double>::min())));
	std::vector<bool> non_null_element (ele_vector_size, true);
	size_t non_null_elements (0);
	for (size_t i(0); i < ele_vector_size; i++) {
		if ((mesh->ele_vector[i])->getVolume() < epsilon) {
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
}

} // end namespace FileIO
