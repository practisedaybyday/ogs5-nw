/*
 * TetGenInterface.cpp
 *
 *  Created on: Sep 12, 2011
 *      Author: TF
 */

#include <cstddef>
#include <string>

// FileIO
#include "MeshIO/TetGenInterface.h"

// Base
#include "StringTools.h"

// MSH
#include "msh_elem.h"
#include "msh_mesh.h"
#include "msh_node.h"

namespace FileIO
{
TetGenInterface::TetGenInterface() :
	_mesh (NULL), _zero_based_idx (true)
{
}

TetGenInterface::~TetGenInterface()
{
}

MeshLib::CFEMesh* TetGenInterface::readTetGenMesh (std::string const& nodes_fname,
                                                   std::string const& ele_fname)
{
	std::ifstream ins_nodes (nodes_fname.c_str());
	std::ifstream ins_ele (ele_fname.c_str());

	if (!ins_nodes || !ins_ele)
	{
		if (!ins_nodes)
			std::cout << "TetGenInterface::readTetGenMesh failed to open " <<
			nodes_fname << std::endl;
		if (!ins_ele)
			std::cout << "TetGenInterface::readTetGenMesh failed to open " <<
			ele_fname << std::endl;
		return NULL;
	}

	_mesh = new MeshLib::CFEMesh();

	if (!readNodesFromStream (ins_nodes))
	{
		delete _mesh;
		return NULL;
	}

	_mesh->InitialNodesNumber();

	if (!readElementsFromStream (ins_ele))
	{
		delete _mesh;
		return NULL;
	}

	return _mesh;
}

bool TetGenInterface::readNodesFromStream (std::ifstream &ins)
{
	std::string line;
	getline (ins, line);
	size_t pos_beg (line.find_first_not_of(" "));
	size_t n_nodes, dim, n_attributes;
	bool boundary_markers;
	bool not_read_header (true);

	while (!ins.fail() && not_read_header)
	{
		line = line.substr(pos_beg);
		if (line.compare(0,1,"#") == 0)
		{
			// this line is a comment - skip
			getline (ins, line);
			pos_beg = line.find_first_not_of(" ");
		}
		else
			// read header line
			not_read_header = !parseNodesFileHeader(line,
			                                        n_nodes,
			                                        dim,
			                                        n_attributes,
			                                        boundary_markers);
	}
	if (not_read_header)
		return false;
	if (!parseNodes(ins, n_nodes, dim))
		return false;

	return true;
}

bool TetGenInterface::parseNodesFileHeader(std::string &line,
                                           size_t& n_nodes,
                                           size_t& dim,
                                           size_t& n_attributes,
                                           bool& boundary_markers) const
{
	size_t pos_beg, pos_end;

	// number of nodes
	pos_beg = line.find_first_not_of (" ");
	pos_end = line.find_first_of(" ", pos_beg);
	if (pos_beg != std::string::npos && pos_end != std::string::npos)
		n_nodes = str2number<size_t> (line.substr(pos_beg, pos_end - pos_beg));
	else
	{
		std::cout <<
		"TetGenInterface::parseNodesFileHeader could not correct read TetGen mesh header - number of nodes"
		          << std::endl;
		return false;
	}
	// dimension
	pos_beg = line.find_first_not_of (" ", pos_end);
	pos_end = line.find_first_of(" ", pos_beg);
	dim = str2number<size_t> (line.substr(pos_beg, pos_end - pos_beg));
	// number of attributes
	pos_beg = line.find_first_not_of (" ", pos_end);
	pos_end = line.find_first_of(" ", pos_beg);
	n_attributes = str2number<size_t> (line.substr(pos_beg, pos_end - pos_beg));
	// boundary marker at nodes?
	pos_beg = line.find_first_not_of (" ", pos_end);
	pos_end = line.find_first_of(" ", pos_beg);
	if (pos_end == std::string::npos)
		pos_end = line.size();
	if ((line.substr(pos_beg, pos_end - pos_beg)).compare("1") == 0)
		boundary_markers = true;
	else
		boundary_markers = false;

	return true;
}

bool TetGenInterface::parseNodes(std::ifstream& ins, size_t n_nodes, size_t dim)
{
	size_t pos_beg, pos_end;
	std::string line;
	double* coordinates (static_cast<double*> (alloca (sizeof(double) * dim)));

	for (size_t k(0); k < n_nodes && !ins.fail(); k++)
	{
		getline (ins, line);
		if (!ins.fail())
		{
			if (!line.empty())
			{
				pos_end = 0;
				// read id
				size_t id;
				pos_beg = line.find_first_not_of(" ", pos_end);
				pos_end = line.find_first_of(" \n", pos_beg);
				if (pos_beg != std::string::npos && pos_end != std::string::npos)
				{
					id =
					        str2number<size_t>(line.substr(pos_beg, pos_end -
					                                       pos_beg));
					if (k == 0 && id == 0)
						_zero_based_idx = true;
				}
				else
				{
					std::cout << "error reading id of node " << k <<
					" in TetGenInterface::parseNodes" << std::endl;
					return false;
				}
				// read coordinates
				for (size_t i(0); i < dim; i++)
				{
					pos_beg = line.find_first_not_of(" ", pos_end);
					pos_end = line.find_first_of(" \n", pos_beg);
					if (pos_end == std::string::npos)
						pos_end = line.size();
					if (pos_beg != std::string::npos)
						coordinates[i] =
						        str2number<double>(line.substr(pos_beg,
						                                       pos_end -
						                                       pos_beg));
					else
					{
						std::cout << "error reading coordinate " << i <<
						" of node " << k <<
						" in TetGenInterface::parseNodes" <<
						std::endl;
						return false;
					}
				}
				if (!_zero_based_idx)
					id--;
				// since CFEMesh is our friend we can access private data of mesh
				_mesh->nod_vector.push_back(new MeshLib::CNode(id, coordinates[0],
				                                               coordinates[1],
				                                               coordinates[2]));
				// read attributes and boundary markers ... - at the moment we do not use this information
			}
		}
		else
		{
			std::cout << "error reading node " << k <<
			" in TetGenInterface::parseNodes" << std::endl;
			return false;
		}
	}

	return true;
}

bool TetGenInterface::readElementsFromStream(std::ifstream &ins)
{
	std::string line;
	getline (ins, line);
	size_t pos_beg (line.find_first_not_of(" "));
	size_t n_tets, n_nodes_per_tet;
	bool region_attributes;
	bool not_read_header (true);

	while (!ins.fail() && not_read_header)
	{
		line = line.substr(pos_beg);
		if (line.compare(0,1,"#") == 0)
		{
			// this line is a comment - skip
			getline (ins, line);
			pos_beg = line.find_first_not_of(" ");
		}
		else
			// read header line
			not_read_header = !parseElementsFileHeader(line,
			                                           n_tets,
			                                           n_nodes_per_tet,
			                                           region_attributes);
	}
	if (not_read_header)
		return false;
	if (!parseElements(ins, n_tets, n_nodes_per_tet, region_attributes))
		return false;

	return true;
}

bool TetGenInterface::parseElementsFileHeader(std::string &line,
                                              size_t& n_tets,
                                              size_t& n_nodes_per_tet,
                                              bool& region_attribute) const
{
	size_t pos_beg, pos_end;

	// number of tetrahedras
	pos_beg = line.find_first_not_of (" ");
	pos_end = line.find_first_of(" ", pos_beg);
	if (pos_beg != std::string::npos && pos_end != std::string::npos)
		n_tets = str2number<size_t> (line.substr(pos_beg, pos_end - pos_beg));
	else
	{
		std::cout <<
		"TetGenInterface::parseElementsFileHeader could not correct read TetGen mesh header - number of tetrahedras"
		          << std::endl;
		return false;
	}
	// nodes per tet - either 4 or 10
	pos_beg = line.find_first_not_of (" \t", pos_end);
	pos_end = line.find_first_of(" \t", pos_beg);
	n_nodes_per_tet = str2number<size_t> (line.substr(pos_beg, pos_end - pos_beg));
	// region attribute at tetrahedra?
	pos_beg = line.find_first_not_of (" \t", pos_end);
	pos_end = line.find_first_of(" \t\n", pos_beg);
	if (pos_end == std::string::npos)
		pos_end = line.size();
	if ((line.substr(pos_beg, pos_end - pos_beg)).compare("1") == 0)
		region_attribute = true;
	else
		region_attribute = false;

	return true;
}

bool TetGenInterface::parseElements(std::ifstream& ins, size_t n_tets, size_t n_nodes_per_tet,
                                    bool region_attribute)
{
	size_t pos_beg, pos_end;
	std::string line;
	size_t* ids (static_cast<size_t*>(alloca (sizeof (size_t) * n_nodes_per_tet)));

	for (size_t k(0); k < n_tets && !ins.fail(); k++)
	{
		getline (ins, line);
		if (!ins.fail())
		{
			if (!line.empty())
			{
				pos_end = 0;
				// read id
				size_t id;
				pos_beg = line.find_first_not_of(" ", pos_end);
				pos_end = line.find_first_of(" \n", pos_beg);
				if (pos_beg != std::string::npos && pos_end != std::string::npos)
					id =
					        str2number<size_t>(line.substr(pos_beg, pos_end -
					                                       pos_beg));
				else
				{
					std::cout << "error reading id of tetrahedra " << k <<
					" in TetGenInterface::parseElements" << std::endl;
					return false;
				}
				// read node ids
				for (size_t i(0); i < n_nodes_per_tet; i++)
				{
					pos_beg = line.find_first_not_of(" ", pos_end);
					pos_end = line.find_first_of(" ", pos_beg);
					if (pos_end == std::string::npos)
						pos_end = line.size();
					if (pos_beg != std::string::npos && pos_end !=
					    std::string::npos)
						ids[i] =
						        str2number<size_t>(line.substr(pos_beg,
						                                       pos_end -
						                                       pos_beg));
					else
					{
						std::cout << "error reading node " << i <<
						" of tetrahedra " << k <<
						" in TetGenInterface::parseElements" <<
						std::endl;
						return false;
					}
				}
				if (!_zero_based_idx)
				{
					id--;
					for (size_t i(0); i < n_nodes_per_tet; i++)
						ids[i]--;
				}
				// since CFEMesh is our friend we can access private data of mesh
				MeshLib::CElem* elem (new MeshLib::CElem(id));
				elem->setElementProperties (MshElemType::TETRAHEDRON, false);
				std::vector<MeshLib::CNode*> ele_nodes(n_nodes_per_tet);
				for (size_t i(0); i < n_nodes_per_tet; i++)
					ele_nodes[i] = _mesh->nod_vector[ids[i]];

				elem->setNodes (ele_nodes);
				_mesh->ele_vector.push_back(elem);
				// read region attribute - this is something like material group
				if (region_attribute)
				{
					pos_beg = line.find_first_not_of(" ", pos_end);
					pos_end = line.find_first_of(" ", pos_beg);
					if (pos_end == std::string::npos)
						pos_end = line.size();
					if (pos_beg != std::string::npos && pos_end !=
					    std::string::npos)
						elem->setPatchIndex (str2number<int>(line.substr(
						                                             pos_beg,
						                                             pos_end
						                                             - pos_beg)));
					else
					{
						std::cout <<
						"error reading region attribute of tetrahedra "
						          << k <<
						" in TetGenInterface::parseElements" <<
						std::endl;
						return false;
					}
				}
			}
		}
		else
		{
			std::cout << "error reading node " << k <<
			" in TetGenInterface::parseElements" << std::endl;
			return false;
		}
	}
	return true;
}
}
