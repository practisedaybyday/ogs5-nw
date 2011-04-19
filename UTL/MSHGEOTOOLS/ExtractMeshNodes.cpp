/*
 * ExtractMeshNodes.cpp
 *
 *  Created on: Dec 3, 2010
 *      Author: TF
 */

#include "ExtractMeshNodes.h"

// BASELIB
#include "quicksort.h"

// GEO
#include "Point.h"

namespace Mesh_Group {

ExtractMeshNodes::ExtractMeshNodes(const CFEMesh* msh) :
	_msh (msh), _gli_pnt_offset (0)
{
}

void ExtractMeshNodes::writeMeshNodeIDs (std::ostream& os, std::ostream& gli_out, const GEOLIB::Polygon& polygon)
{
	// get all nodes of mesh
	const std::vector<Mesh_Group::CNode*>& msh_nodes (_msh->getNodeVector());

	std::vector<size_t> node_indices;

	for (size_t j(0); j<msh_nodes.size(); j++) {
		if (msh_nodes[j]->Interior()) {
			GEOLIB::Point pnt (msh_nodes[j]->X(), msh_nodes[j]->Y(), 0);
			if (polygon.isPntInPolygon(pnt)) {
				node_indices.push_back (j);
			}
		}
	}
	// write data
	for (size_t k(0); k<node_indices.size(); k++) {
		os << node_indices[k] << std::endl;
	}

	for (size_t k(0); k<node_indices.size(); k++) {
		gli_out << k + _gli_pnt_offset << " " << msh_nodes[node_indices[k]]->X() << " " <<  msh_nodes[node_indices[k]]->Y() << " " << msh_nodes[node_indices[k]]->Z() << std::endl;
	}
	_gli_pnt_offset += node_indices.size();
}

void ExtractMeshNodes::writeTopSurfaceMeshNodeIDs (std::ostream& os, std::ostream& gli_out, const GEOLIB::Polygon& polygon)
{
	// get all nodes of mesh
	const std::vector<Mesh_Group::CNode*>& msh_nodes (_msh->getNodeVector());

	std::vector<GEOLIB::PointWithID> nodes_as_points;

	for (size_t j(0); j<msh_nodes.size(); j++) {
//		if (msh_nodes[j]->Interior()) {
			GEOLIB::Point pnt (msh_nodes[j]->X(), msh_nodes[j]->Y(), 0);
			if (polygon.isPntInPolygon(pnt)) {
				nodes_as_points.push_back (GEOLIB::PointWithID (msh_nodes[j]->X(), msh_nodes[j]->Y(), msh_nodes[j]->Z(), j));
			}
//		}
	}

	std::vector<size_t> perm;
	for (size_t k(0); k<nodes_as_points.size(); k++) {
		perm.push_back(k);
	}
	Quicksort<GEOLIB::PointWithID> (nodes_as_points, 0, nodes_as_points.size(), perm);

	double eps (sqrt(std::numeric_limits<double>::min()));

	// write data
	for (size_t k(1); k<nodes_as_points.size(); k++) {
		const GEOLIB::PointWithID& p0 (nodes_as_points[k-1]);
		const GEOLIB::PointWithID& p1 (nodes_as_points[k]);
		if (fabs (p0[0]-p1[0]) > eps || fabs (p0[1]-p1[1]) > eps) {
			os << p0.getID() << std::endl;
		}
	}
	// write last point
	os << nodes_as_points[nodes_as_points.size()-1].getID() << std::endl;

	size_t n_nodes (0);
	gli_out.precision (14);
	for (size_t k(1); k<nodes_as_points.size(); k++) {
		const GEOLIB::PointWithID& p0 (nodes_as_points[k-1]);
		const GEOLIB::PointWithID& p1 (nodes_as_points[k]);
		if (fabs (p0[0]-p1[0]) > eps || fabs (p0[1]-p1[1]) > eps) {
			gli_out << n_nodes + _gli_pnt_offset << " " << std::scientific << p0 << " $NAME " << p0.getID() << std::endl;
			n_nodes++;
		}
	}
	// write last point
	gli_out << n_nodes + _gli_pnt_offset << " " << std::scientific << nodes_as_points[nodes_as_points.size()-1] << " $NAME " << nodes_as_points[nodes_as_points.size()-1].getID() << std::endl;
	n_nodes++;
	_gli_pnt_offset += n_nodes;
}

void ExtractMeshNodes::writeMesh2DNodeIDAndArea (std::ostream& os, std::ostream& gli_out, const GEOLIB::Polygon& polygon)
{
	// get all nodes of mesh
	const std::vector<Mesh_Group::CNode*>& msh_nodes (_msh->getNodeVector());

	// store node id
	std::vector<size_t> node_ids;

	for (size_t j(0); j<msh_nodes.size(); j++) {
		GEOLIB::Point pnt (msh_nodes[j]->X(), msh_nodes[j]->Y(), 0);
		if (polygon.isPntInPolygon(pnt)) {
			node_ids.push_back (j);
		}
	}
	std::sort (node_ids.begin(), node_ids.end());

	size_t n_nodes (node_ids.size());
	std::vector<double> areas (n_nodes, 0.0);
	// in order to compute the area we need the mesh elements
	const std::vector<Mesh_Group::CElem*>& msh_elem (_msh->getElementVector());
	for (size_t k(0); k<n_nodes; k++) {
		// get all associated mesh elements
		std::vector<size_t> const & mesh_elem_ids(msh_nodes[node_ids[k]]->getConnectedElementIDs());
		size_t n_mesh_elem_ids (mesh_elem_ids.size());
		// get areas for mesh elements
		double area (0.0);
		for (size_t j(0); j<n_mesh_elem_ids; j++) {
			// check if all mesh nodes of the element are inside the polygon
			std::vector<size_t> mesh_element_node_indices;
			msh_elem[mesh_elem_ids[j]]->getNodeIndices (mesh_element_node_indices);
			size_t n_of_nodes_of_element (mesh_element_node_indices.size());
			bool found (true);
			for (size_t i(0); i<n_of_nodes_of_element && found; i++) {
				std::vector<size_t>::iterator it(std::find (node_ids.begin(), node_ids.end(), mesh_element_node_indices[i]));
				if (it == node_ids.end())
					found = false;
			}
			if (found)
				area += msh_elem[mesh_elem_ids[j]]->calcVolume();
		}
		areas[k] = area/3.0;
	}

	// write ids and areas
	for (size_t k(0); k<n_nodes; k++) {
		if (areas[k] > std::numeric_limits<double>::epsilon())
			os << node_ids[k] << " " << areas[k] << std::endl;
	}

	n_nodes = 0;
	gli_out.precision (14);
	for (std::vector<size_t>::const_iterator it (node_ids.begin()); it != node_ids.end(); it++) {
		gli_out << n_nodes + _gli_pnt_offset << " " << std::scientific
			<< msh_nodes[*it]->X() << " " << msh_nodes[*it]->Y()<< " " << msh_nodes[*it]->Z()
			<< " $NAME " << *it << std::endl;
		n_nodes++;
	}
	_gli_pnt_offset += n_nodes;
}

void ExtractMeshNodes::writeNearestMeshNodeToPoint (std::ostream& os, std::ostream& gli_out, GEOLIB::Point const & pnt)
{
	size_t node_id (_msh->GetNODOnPNT (&pnt));
	os << node_id << std::endl;

	gli_out.precision (14);
	gli_out << _gli_pnt_offset << " " << std::scientific
		<< (_msh->getNodeVector()[node_id])->X() << " "
		<< (_msh->getNodeVector()[node_id])->Y() << " "
		<< (_msh->getNodeVector()[node_id])->Z() << std::endl;

	_gli_pnt_offset++;
}

} // end namespace Mesh_Group
