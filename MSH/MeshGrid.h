/*
 * MeshGrid.h
 *
 *  Created on: Feb 2, 2012
 *      Author: TF
 */

#ifndef MESHGRID_H_
#define MESHGRID_H_

#include <vector>

// GeoLib
#include "AxisAlignedBoundingBox.h"

// MeshLib
#include "msh_mesh.h"
#include "msh_node.h"

namespace MeshLib {

class MeshGrid : public GEOLIB::AABB
{
public:
	MeshGrid(MeshLib::CFEMesh const& mesh, size_t max_num_per_grid_cell = 500);
	virtual ~MeshGrid();

	size_t getIndexOfNearestNode(double const*const pnt) const;

#ifndef NDEBUG
	/**
	 * Method creates a geometry for every mesh grid box. Additionally it
	 * creates one geometry containing all the box geometries.
	 * @param geo_obj
	 */
	void createMeshGridGeometry(GEOLIB::GEOObjects* geo_obj) const;
#endif

private:
	inline void getGridCoords(double const*const pnt, size_t* coords) const;
	bool calcNearestNodeInGrid(double const* const pnt, size_t const* const coords,
					double &sqr_min_dist, size_t &global_idx) const;
	double _step_sizes[3];
	double _inverse_step_sizes[3];
	size_t _n_steps[3];
	std::vector<MeshLib::CNode*>* _grid_quad_to_node_map;
};

} // end namespace MeshLib

#endif /* MESHGRID_H_ */
