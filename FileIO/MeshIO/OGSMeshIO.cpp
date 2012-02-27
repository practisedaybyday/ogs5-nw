/*
 * OGSMeshIO.cpp
 *
 *  Created on: Mar 3, 2011
 */

#include "GEOObjects.h"
#include "MeshIO/OGSMeshIO.h"
#include "msh_lib.h"
#include "msh_mesh.h"
#include "msh_elem.h"

namespace FileIO
{

OGSMeshIO::OGSMeshIO()
: _mesh(NULL)
{
}

MeshLib::CFEMesh* OGSMeshIO::loadMeshFromFile(std::string const& fname)
{
	std::cout << "FEMRead ... " << std::flush;
/*
   #ifndef NDEBUG
    QTime myTimer;
    myTimer.start();
   #endif
 */
	std::vector<MeshLib::CFEMesh*> mesh_vec;
	FEMRead(fname.substr(0, fname.length() - 4), mesh_vec);
	if (!mesh_vec.empty())
	{
/*
   #ifndef NDEBUG
        QTime constructTimer;
        constructTimer.start();
   #endif
 */
		mesh_vec[mesh_vec.size() - 1]->ConstructGrid();

		std::cout << "Nr. Nodes: " <<
		mesh_vec[mesh_vec.size() - 1]->nod_vector.size() << std::endl;
		std::cout << "Nr. Elements: " <<
		mesh_vec[mesh_vec.size() - 1]->ele_vector.size() << std::endl;
/*
   #ifndef NDEBUG
        std::cout << "constructGrid time: " << constructTimer.elapsed() << " ms" << std::endl;
   #endif
 */
		mesh_vec[mesh_vec.size() - 1]->FillTransformMatrix();
/*
   #ifndef NDEBUG
        std::cout << "Loading time: " << myTimer.elapsed() << " ms" << std::endl;
   #endif
 */
		return mesh_vec[mesh_vec.size() - 1];
	}

	std::cout << "Failed to load the mesh file: " << fname << std::endl;
	return NULL;
}

void OGSMeshIO::write(std::ostream &out)
{
	if(!_mesh)
	{
		std::cout << "OGSMeshIO cannot write: no mesh set!" << std::endl;
		return;
	}
	
	out << "#FEM_MSH" << std::endl;

	out << "$PCS_TYPE" << std::endl << "  " << _mesh->pcs_name << std::endl;

	out << "$NODES" << std::endl << "  ";
	const size_t n_nodes(_mesh->GetNodesNumber(false));
	out << n_nodes << std::endl;
	for (size_t i(0); i < n_nodes; i++)
	{
		double const* const coords (_mesh->nod_vector[i]->getData());
		out << i << " " << coords[0] << " " << coords[1] << " " << coords[2] << std::endl;
	}

	out << "$ELEMENTS" << std::endl << "  ";
	writeElementsExceptLines(_mesh->ele_vector, out);
//	const size_t ele_vector_size (_mesh->ele_vector.size());
//	const double epsilon (std::numeric_limits<double>::epsilon());
//	std::vector<bool> non_null_element (ele_vector_size, true);
//	size_t non_null_elements (0);
//	for (size_t i(0); i < ele_vector_size; i++)
//	{
//		if ((_mesh->ele_vector[i])->calcVolume() < epsilon) // || (_mesh->ele_vector[i])->GetElementType() == MshElemType::LINE) {
//			non_null_element[i] = false;
//		else
//			non_null_elements++;
//	}
//	out << non_null_elements << std::endl;
//	for (size_t i(0), k(0); i < ele_vector_size; i++) {
//		if (non_null_element[i])
//		{
//			_mesh->ele_vector[i]->SetIndex(k);
//			_mesh->ele_vector[i]->WriteIndex(out);
//			k++;
//		}
//	}

	out << " $LAYER" << std::endl;
	out << "  ";
	out << _mesh->_n_msh_layer << std::endl;
	out << "#STOP" << std::endl;
}

void OGSMeshIO::writeElementsExceptLines(std::vector<MeshLib::CElem*> const& ele_vec, std::ostream &out)
{
	const size_t ele_vector_size (ele_vec.size());
	std::vector<bool> non_line_element (ele_vector_size, true);
	size_t n_non_line_elements(0);
	for (size_t i(0); i < ele_vector_size; i++) {
		if ((ele_vec[i])->GetElementType() == MshElemType::LINE) {
			non_line_element[i] = false;
		} else {
			n_non_line_elements++;
		}
	}
	out << n_non_line_elements << std::endl;
	for (size_t i(0), k(0); i < ele_vector_size; i++) {
		if (non_line_element[i]) {
			const size_t tmp_idx (ele_vec[i]->GetIndex());
			(const_cast<MeshLib::CElem*>(ele_vec[i]))->SetIndex(k);
			ele_vec[i]->WriteIndex(out);
			(const_cast<MeshLib::CElem*>(ele_vec[i]))->SetIndex(tmp_idx);
			k++;
		}
	}
}

void OGSMeshIO::setMesh(MeshLib::CFEMesh const* mesh)
{
	_mesh = mesh;
}

} // end namespace FileIO
