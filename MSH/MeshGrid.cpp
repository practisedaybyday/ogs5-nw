/*
 * MeshGrid.cpp
 *
 *  Created on: Feb 2, 2012
 *      Author: TF
 */

#include "MeshGrid.h"

namespace MeshLib {

MeshGrid::MeshGrid(MeshLib::CFEMesh const& mesh) :
	GEOLIB::AABB(), _grid_quad_to_node_map(NULL)
{
	// compute axis aligned bounding box
	std::vector<MeshLib::CNode*> const& nodes(mesh.getNodeVector());
	const size_t n_nodes(mesh.GetNodesNumber(false));
	for (size_t k(0); k<n_nodes; k++) {
		this->update(nodes[k]->getData());
	}

	for (size_t k(0); k<3; k++) {
		// make the bounding box a little bit bigger,
		// such that the node with maximal coordinates fits into the grid
		_max_pnt[k] *= (1.0+1e-6);
		if (fabs(_max_pnt[k]) < std::numeric_limits<double>::epsilon()) {
			_max_pnt[k] = (_max_pnt[k] - _min_pnt[k]) * (1.0+1e-6);
		}
		_delta[k] = _max_pnt[k] - _min_pnt[k];
	}

	// *** condition: n_nodes / (_n_steps[0] * _n_steps[1] * _n_steps[2]) < 500
	// *** with _n_steps[1] = _n_steps[0] * _delta[1]/_delta[0], _n_steps[2] = _n_steps[0] * _delta[2]/_delta[0]
	if (fabs(_delta[1]) < std::numeric_limits<double>::epsilon() || fabs(_delta[2]) < std::numeric_limits<double>::epsilon()) {
		// 1d case y = z = 0
		if (fabs(_delta[1]) < std::numeric_limits<double>::epsilon() && fabs(_delta[2]) < std::numeric_limits<double>::epsilon()) {
			_n_steps[0] = static_cast<size_t>(ceil(n_nodes / 500.0));
			_n_steps[1] = 1;
			_n_steps[2] = 1;
		} else {
			// 1d case x = z = 0
			if (fabs(_delta[0]) < std::numeric_limits<double>::epsilon() && fabs(_delta[2]) < std::numeric_limits<double>::epsilon()) {
				_n_steps[0] = 1;
				_n_steps[1] = static_cast<size_t>(ceil(n_nodes / 500.0));
				_n_steps[2] = 1;
			} else {
				// 1d case x = y = 0
				if (fabs(_delta[0]) < std::numeric_limits<double>::epsilon() && fabs(_delta[1]) < std::numeric_limits<double>::epsilon()) {
					_n_steps[0] = 1;
					_n_steps[1] = 1;
					_n_steps[2] = static_cast<size_t>(ceil(n_nodes / 500.0));
				} else {
					// 2d case
					if (fabs(_delta[1]) < std::numeric_limits<double>::epsilon()) {
						_n_steps[0] = static_cast<size_t>(ceil(sqrt(n_nodes * _delta[0] / (500*_delta[2]))));
						_n_steps[1] = 1;
						_n_steps[2] = static_cast<size_t>(ceil(_n_steps[0] * _delta[2] / _delta[0]));
					} else {
						_n_steps[0] = static_cast<size_t>(ceil(sqrt(n_nodes * _delta[0] / (500*_delta[1]))));
						_n_steps[1] = static_cast<size_t>(ceil(_n_steps[0] * _delta[1] / _delta[0]));
						_n_steps[2] = 1;
					}
				}
			}
		}
	} else {
		// 3d case
		_n_steps[0] = static_cast<size_t>(ceil(pow(n_nodes * _delta[0]*_delta[0] / (500*_delta[1]*_delta[2]), 1. / 3.)));
		_n_steps[1] = static_cast<size_t>(ceil(_n_steps[0] * _delta[1] / _delta[0]));
		_n_steps[2] = static_cast<size_t>(ceil(_n_steps[0] * _delta[2] / _delta[0]));
	}

	const size_t n_plane (_n_steps[0]*_n_steps[1]);
	_grid_quad_to_node_map = new std::vector<MeshLib::CNode*> [n_plane*_n_steps[2]];

	// some frequently used expressions to fill the grid vectors
	for (size_t k(0); k<3; k++) {
		_step_sizes[k] = _delta[k] / _n_steps[k];
		_inverse_step_sizes[k] = 1.0 / _step_sizes[k];
	}

	// fill the grid vectors
	for (size_t l(0); l<n_nodes; l++) {
		double const*const node(nodes[l]->getData());
		const size_t i ((node[0]-_min_pnt[0]) * _inverse_step_sizes[0]);
		const size_t j ((node[1]-_min_pnt[1]) * _inverse_step_sizes[1]);
		const size_t k ((node[2]-_min_pnt[2]) * _inverse_step_sizes[2]);

		if (i >= _n_steps[0] || j >= _n_steps[1] || k >= _n_steps[2]) {
			std::cout << "error computing indices " << std::endl;
		}

		_grid_quad_to_node_map[i + j*_n_steps[0]+k*n_plane].push_back (nodes[l]);
	}

#ifndef NDEBUG
	size_t nodes_cnt(0);
	for (size_t k(0); k<n_plane*_n_steps[2]; k++)
		nodes_cnt += _grid_quad_to_node_map[k].size();

	std::cout << "mesh has " << n_nodes << " nodes, grid has " << nodes_cnt << " nodes" << std::endl;
	assert(n_nodes==nodes_cnt);
#endif
}

MeshGrid::~MeshGrid()
{
	delete [] _grid_quad_to_node_map;
}

void MeshGrid::getGridCornerPoints(double const*const node, double* llf, double* urb) const
{
	size_t coords[3];
	getGridCoords(node, coords);
	for (unsigned l(0); l<3; l++) {
		llf[l] = _min_pnt[l] + coords[l] * _step_sizes[l];
		urb[l] = _min_pnt[l] + (coords[l]+1) * _step_sizes[l];
	}
}

std::vector<MeshLib::CNode*> const& MeshGrid::getNodesInGrid(double const*const node) const
{
	size_t coords[3];
	getGridCoords(node, coords);
	return _grid_quad_to_node_map[coords[0] + coords[1]*_n_steps[0]+coords[2]*_n_steps[0]*_n_steps[1]];
}

std::vector<MeshLib::CNode*> const& MeshGrid::getNodesInGrid(size_t const*const coords) const
{
	return _grid_quad_to_node_map[coords[0] + coords[1]*_n_steps[0]+coords[2]*_n_steps[0]*_n_steps[1]];
}

void MeshGrid::getGridCoords(double const*const node, size_t* coords) const
{
	for (size_t k(0); k<3; k++)
		coords[k] = (node[k]-_min_pnt[k]) * _inverse_step_sizes[k];
}

size_t MeshGrid::getIndexOfNearestNode(double const*const pnt) const
{
	size_t coords[3];
	getGridCoords(pnt, coords);

	double sqr_min_dist (MathLib::sqrDist(&_min_pnt, &_max_pnt));
	size_t global_idx(std::numeric_limits<size_t>::max());
	if (calcNearestNodeInGrid(pnt, coords, sqr_min_dist, global_idx)) {
		double sqr_min_dist_tmp;
		size_t global_idx_tmp;

		// check all border cuboids
		size_t tmp_coords[3];
		for (size_t i(0); i<3; i++) {
			tmp_coords[0] = coords[0]-1+i;
			for (size_t j(0); j<3; j++) {
				tmp_coords[1] = coords[1]-1+j;
				for (size_t k(0); k<3; k++) {
					tmp_coords[2] = coords[2]-1+k;
					if (!(i == j && i == k)) {
						if (calcNearestNodeInGrid(pnt, tmp_coords, sqr_min_dist_tmp, global_idx_tmp)) {
							if (sqr_min_dist_tmp < sqr_min_dist) {
								sqr_min_dist = sqr_min_dist_tmp;
								global_idx = global_idx_tmp;
							}
						}
					}
				}  // end k
			} // end j
		} // end i
	}

	return global_idx;
}

bool MeshGrid::calcNearestNodeInGrid(double const* const pnt, size_t const* const coords,
				double &sqr_min_dist, size_t &global_idx) const
{
	if (coords[0] >= _n_steps[0] || coords[1] >= _n_steps[1] || coords[2] >= _n_steps[2])
		return false;

	const size_t grid_idx (coords[0] + coords[1] * _n_steps[0] + coords[2] * _n_steps[0] * _n_steps[1]);
	std::vector<MeshLib::CNode*> const& nodes(_grid_quad_to_node_map[grid_idx]);
	if (nodes.empty()) return false;

	const size_t n_nodes(nodes.size());
	sqr_min_dist = MathLib::sqrDist(nodes[0]->getData(), pnt);
	global_idx = nodes[0]->GetIndex();
	for (size_t i(1); i < n_nodes; i++) {
		const double sqr_dist(MathLib::sqrDist(nodes[i]->getData(), pnt));
		if (sqr_dist < sqr_min_dist) {
			sqr_min_dist = sqr_dist;
			global_idx = nodes[i]->GetIndex();
		}
	}
	return true;
}

void MeshGrid::getNodeVectorsInAxisAlignedBoundingBox(GEOLIB::Point const& ll,
				GEOLIB::Point const& ur, size_t &n_node_vectors, std::vector<MeshLib::CNode*> * & node_vectors)
{
	n_node_vectors = _n_steps[0]*_n_steps[1]*_n_steps[2];
	node_vectors = new std::vector<MeshLib::CNode*>[n_node_vectors];
	for (size_t k(0); k<n_node_vectors; k++)
		node_vectors[k] = _grid_quad_to_node_map[k];
}

} // end namespace MeshLib
