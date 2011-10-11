/**
 * \file MshEditor.cpp
 * 2011/06/15 KR Initial implementation
 */

#include "MshEditor.h"
#include "msh_mesh.h"

MeshLib::CFEMesh* MshEditor::removeMeshNodes(MeshLib::CFEMesh* mesh,
                                             const std::vector<size_t> &nodes)
{
	MeshLib::CFEMesh* new_mesh (new MeshLib::CFEMesh(*mesh));

	// delete nodes and their connected elements and replace them with null pointers
	size_t delNodes = nodes.size();
	for (size_t i = 0; i < delNodes; i++)
	{
		MeshLib::CNode* node = new_mesh->nod_vector[nodes[i]];
		std::vector<size_t> conn_elems = node->getConnectedElementIDs();
		for (size_t j = 0; j < conn_elems.size(); j++)
		{
			delete new_mesh->ele_vector[conn_elems[j]];
			new_mesh->ele_vector[conn_elems[j]] = NULL;
		}
		delete new_mesh->nod_vector[nodes[i]];
		new_mesh->nod_vector[nodes[i]] = NULL;
	}

	// create map to adjust node indices in element vector
	size_t nNodes = new_mesh->nod_vector.size();
	std::vector<int> id_map;
	size_t count = 0;
	for (size_t i = 0; i < nNodes; i++)
	{
		if (new_mesh->nod_vector[i])
		{
			new_mesh->nod_vector[i]->SetIndex(count);
			id_map.push_back(count);
			count++;
		}
		else
			id_map.push_back(-1);
	}

	// erase null pointers from node- and element vectors
	for (std::vector<MeshLib::CElem*>::iterator it = new_mesh->ele_vector.begin();
	     it != new_mesh->ele_vector.end(); )
	{
		if (*it)
			++it;
		else
			it = new_mesh->ele_vector.erase(it);
	}

	for (std::vector<MeshLib::CNode*>::iterator it = new_mesh->nod_vector.begin();
	     it != new_mesh->nod_vector.end(); )
	{
		if (*it)
			++it;
		else
			it = new_mesh->nod_vector.erase(it);
	}

	// re-adjust node indices
	size_t nElems = new_mesh->ele_vector.size();
	for (size_t i = 0; i < nElems; i++)
	{
		MeshLib::CElem* elem = new_mesh->ele_vector[i];
		size_t nElemNodes = elem->GetNodesNumber(false);
		for (size_t j = 0; j < nElemNodes; j++)
			elem->SetNodeIndex(j, id_map[elem->GetNodeIndex(j)]);
	}

	return new_mesh;
}

