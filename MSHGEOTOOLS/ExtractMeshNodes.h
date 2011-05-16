/*
 * ExtractMeshNodes.h
 *
 *  Created on: Dec 3, 2010
 *      Author: TF
 */

#ifndef EXTRACTMESHNODES_H_
#define EXTRACTMESHNODES_H_

#include <iostream>

// MSH
#include "msh_mesh.h"
#include "msh_lib.h"

// GEO
#include "GEOObjects.h"
#include "Polygon.h"
#include "PointWithID.h"

namespace Mesh_Group {

/**
 * This class implements an algorithm to extract mesh node ids from a given (extruded) mesh.
 */
class ExtractMeshNodes {
public:
	/**
	 * constructor - take the mesh
	 * @param msh an instance of class CFEMesh
	 */
	ExtractMeshNodes(const CFEMesh* msh);
	/**
	 * This method first projects a mesh node into the x-y-plane (z=0).
	 * Then it checks if this mesh node is within the given polygon
	 * (polygon have to be located in the x-y-plane (z=0).
	 * The IDs of all (projected) mesh nodes within the polygon will
	 * be written to the stream os. For visual control the mesh nodes
	 * will be written (as gli points) to the stream gli_out.
	 * @param os output stream for IDs
	 * @param gli_out output stream for points
	 * @param polygon the polygon that have to be located in the x-y-plane (z=0)
	 */
	void writeMeshNodeIDs (std::ostream& os, std::ostream& gli_out, const GEOLIB::Polygon& polygon);
	/**
	 * This method first projects all mesh nodes into the x-y-plane (z=0).
	 * Then it checks if mesh nodes are within the given polygon
	 * (polygon have to be located in the x-y-plane (z=0). All the mesh
	 * nodes are located within a "cylinder". The method sorts the mesh nodes
	 * lexicographical (first by x, then by y and in the end by z). The id of the
	 * mesh node with largest z coordinate and identical x and y coordinates
	 * will be written to the stream os. For visual control the associated mesh
	 * node will be written to the stream gli_out (as gli point).
	 * @param os output stream for IDs
	 * @param gli_out output stream for points
	 * @param polygon the polygon that have to be located in the x-y-plane (z=0)
	 */
	void writeTopSurfaceMeshNodeIDs (std::ostream& os, std::ostream& gli_out, const GEOLIB::Polygon& polygon);

	void writeMesh2DNodeIDAndArea (std::ostream& os, std::ostream& gli_out, const GEOLIB::Polygon& polygon);

	void writeNearestMeshNodeToPoint (std::ostream& os, std::ostream& gli_out, GEOLIB::Point const & pnt);

	/**
	 * computes the bottom mesh nodes along a polyline
	 * @param ply along the polyline ply
	 * @param bottom_points the bottom mesh nodes as points
	 */
	void getBottomMeshNodesAlongPolylineAsPoints (const GEOLIB::Polyline& ply, std::vector<GEOLIB::Point*>& bottom_points) const;

	/**
	 * Method computes the polygon to a given polyline that is consisting of this polyline, a
	 * projection of this polyline to the bottom of the mesh and the links between these two polylines.
	 * @param polyline the ("defining") polyline which is a part of the polygon,
	 * 	the other part will be created by projection of the polyline to the bottom of the mesh.
	 * @param point_vec point vector that holds (at least) the points of the polygon
	 * @param polygon pointer to the resulting polygon
	 * 	warning: the pointer to an already existing polygon will be destroyed
	 */
	void getPolygonFromPolyline (const GEOLIB::Polyline& polyline, std::vector<GEOLIB::Point*>& point_vec, GEOLIB::Polygon* &polygon) const;

private:
	const CFEMesh* _msh;
	/**
	 * offset for gli point index
	 */
	size_t _gli_pnt_offset;
};

}

#endif /* EXTRACTMESHNODES_H_ */
