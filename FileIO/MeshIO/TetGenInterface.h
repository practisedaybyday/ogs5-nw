/*
 * TetGenInterface.h
 *
 *  Created on: Sep 12, 2011
 *      Author: TF
 */

#ifndef TETGENINTERFACE_H_
#define TETGENINTERFACE_H_

// forward declaration of mesh class
namespace MeshLib {
class CFEMesh;
}

namespace FileIO {

class TetGenInterface {
public:
	TetGenInterface();
	virtual ~TetGenInterface();

	MeshLib::CFEMesh* readTetGenMesh (std::string const& nodes_fname, std::string const& ele_fname);
	/** in order to have a direct access to the
	 * data structures for nodes and elements we make
	 * class TetGenInterface a friend of the mesh class
	 */
	friend class MeshLib::CFEMesh;

private:
	/**
	 * Method reads the nodes from stream and stores them in the node vector of the mesh class.
	 * @param input the input stream
	 */
	bool readNodesFromStream(std::ifstream &input);
	bool parseNodesFileHeader(std::string &line, size_t& n_nodes, size_t& dim,
			size_t& n_attributes, bool& boundary_markers) const;
	bool parseNodes(std::ifstream& ins, size_t n_nodes, size_t dim,
			size_t n_attributes, bool boundary_markers);

	bool readElementsFromStream(std::ifstream &input);
	bool parseElementsFileHeader(std::string &line, size_t& n_tets, size_t& n_nodes_per_tet,
			bool& region_attribute) const;
	bool parseElements(std::ifstream& ins, size_t n_tets, size_t n_nodes_per_tet,
			bool region_attribute);


	MeshLib::CFEMesh* _mesh;
	bool _zero_based_idx;
};

}

#endif /* TETGENINTERFACE_H_ */
