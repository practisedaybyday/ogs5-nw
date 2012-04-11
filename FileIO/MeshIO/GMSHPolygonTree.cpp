/*
 * GMSHPolygonTree.cpp
 *
 *  Created on: Mar 27, 2012
 *      Author: fischeth
 */

#include "MeshIO/GMSHPolygonTree.h"

#include "GMSHNoMeshDensity.h"
#include "GMSHFixedMeshDensity.h"
#include "GMSHAdaptiveMeshDensity.h"

namespace FileIO {

GMSHPolygonTree::GMSHPolygonTree(GEOLIB::Polygon* polygon, GMSHPolygonTree* parent,
				GEOLIB::GEOObjects &geo_objs, std::string const& geo_name,
				GMSHMeshDensityStrategy * mesh_density_strategy) :
	GEOLIB::SimplePolygonTree(polygon, parent), _geo_objs(geo_objs), _geo_name(geo_name),
	_mesh_density_strategy(mesh_density_strategy)
{}

GMSHPolygonTree::~GMSHPolygonTree()
{}

bool GMSHPolygonTree::insertStation(GEOLIB::Point const* station)
{

	if (_node_polygon->isPntInPolygon(*station)) {
		// try to insert station into the child nodes
		for (std::list<SimplePolygonTree*>::const_iterator it (_childs.begin());
			 it != _childs.end(); it++) {
			if (((*it)->getPolygon())->isPntInPolygon (*station)) {
				return dynamic_cast<GMSHPolygonTree*>((*it))->insertStation (station);
			}
		}
		// station did not fit into child nodes -> insert the station into this node
		_stations.push_back (station);
		return true;
	} else {
		return false;
	}
}

void GMSHPolygonTree::insertPolyline (GEOLIB::PolylineWithSegmentMarker * ply)
{
	if (_node_polygon->isPartOfPolylineInPolygon(*ply)) {
		// check childs
		for (std::list<SimplePolygonTree*>::const_iterator it (_childs.begin());
			it != _childs.end(); it++) {
			dynamic_cast<GMSHPolygonTree*>((*it))->insertPolyline (ply);
		}
		_plys.push_back(ply);

		// calculate possible intersection points
		const size_t n_ply_pnts(ply->getNumberOfPoints()-1);
		for (size_t k(0); k<n_ply_pnts; k++) {
			if (! ply->isSegmentMarked(k)) {
				size_t seg_num(0);
				GEOLIB::Point* s(_node_polygon->getIntersectionPointPolygonLine(*(ply->getPoint(k)), *(ply->getPoint(k+1)), seg_num));
				if (s) {
					const double eps (std::numeric_limits<float>::epsilon());

					// *** begin test code
					double max_x_polygon_segment ((*(_node_polygon->getPoint(seg_num)))[0]);
					double min_x_polygon_segment ((*(_node_polygon->getPoint(seg_num+1)))[0]);
					double max_y_polygon_segment ((*(_node_polygon->getPoint(seg_num)))[1]);
					double min_y_polygon_segment ((*(_node_polygon->getPoint(seg_num+1)))[1]);

					double max_x_polyline_segment ((*(ply->getPoint(k)))[0]);
					double min_x_polyline_segment ((*(ply->getPoint(k+1)))[0]);
					double max_y_polyline_segment ((*(ply->getPoint(k)))[1]);
					double min_y_polyline_segment ((*(ply->getPoint(k+1)))[1]);

					if (max_x_polygon_segment > min_x_polygon_segment) {
						std::swap (max_x_polygon_segment, min_x_polygon_segment);
					}
					if (max_y_polygon_segment > min_y_polygon_segment) {
						std::swap (max_y_polygon_segment, min_y_polygon_segment);
					}

					if (max_x_polyline_segment > min_x_polyline_segment) {
						std::swap (max_x_polyline_segment, min_x_polyline_segment);
					}
					if (max_y_polyline_segment > min_y_polyline_segment) {
						std::swap (max_y_polyline_segment, min_y_polyline_segment);
					}

					if ((*s)[0] < min_x_polygon_segment || (*s)[0] > max_x_polygon_segment || (*s)[1] < min_y_polygon_segment || (*s)[1] > max_y_polygon_segment ||
									(*s)[0] < min_x_polyline_segment || (*s)[0] > max_x_polyline_segment || (*s)[1] < min_y_polyline_segment || (*s)[1] > max_y_polyline_segment) {
						std::cout << "[" << *(_node_polygon->getPoint(seg_num)) << "; " << *(_node_polygon->getPoint(seg_num+1)) << "] - [" << *(ply->getPoint(k)) << "; " << *(ply->getPoint(k+1)) << "]" << " -> s: " << *s << std::endl;

						GEOLIB::Point *s1(_node_polygon->getIntersectionPointPolygonLine(*(ply->getPoint(k)), *(ply->getPoint(k+1)), seg_num));
						delete s1;
					}
					// *** end test code
					if (MathLib::sqrDist(s, ply->getPoint(k)) < eps * _node_polygon->getLength(0)) {
						_node_polygon->insertPoint(seg_num+1, ply->getPointID(k));
						ply->markSegment(k, true);
					} else {
						if (MathLib::sqrDist(s, ply->getPoint(k+1)) < eps * _node_polygon->getLength(0)) {
							_node_polygon->insertPoint(seg_num+1, ply->getPointID(k+1));
							ply->markSegment(k, true);
						} else {
							// insert intersection point to point vector of GEOObjects instance
							size_t pnt_id;
							_geo_objs.appendPoint(s, _geo_name, pnt_id);
							// modify the polygon
							_node_polygon->insertPoint(seg_num+1, pnt_id);
							// modify the polyline
							ply->insertPoint(k+1, pnt_id);
							// mark the segment
							if (_node_polygon->isPntInPolygon(*(ply->getPoint(k)))) {
								ply->markSegment(k, true);
							} else {
								ply->markSegment(k+1, true);
							}
							// insert line segment as constraint
							_gmsh_lines_for_constraints.push_back(new GMSHLine(ply->getPointID(k), ply->getPointID(k+1)));
						}
					}
				} else {
					if (_node_polygon->isPntInPolygon(*(ply->getPoint(k)))) {
						ply->markSegment(k, true);
						// line segment is inside polygon
						_gmsh_lines_for_constraints.push_back(new GMSHLine(ply->getPointID(k), ply->getPointID(k+1)));
					}
				}
			}
		}
	}
}

void GMSHPolygonTree::initMeshDensityStrategy()
{
	if (dynamic_cast<GMSHAdaptiveMeshDensity*> (_mesh_density_strategy)) {
		// collect points
		std::vector<GEOLIB::Point const*> pnts;
		const size_t n_pnts_polygon (_node_polygon->getNumberOfPoints());
		for (size_t k(0); k<n_pnts_polygon; k++) {
			pnts.push_back(_node_polygon->getPoint(k));
		}
		getPointsFromSubPolygons(pnts);
		// give collected points to the mesh density strategy
		_mesh_density_strategy->init(pnts);
		// insert constraints
		dynamic_cast<GMSHAdaptiveMeshDensity*>(_mesh_density_strategy)->addPoints(_stations);
		std::vector<GEOLIB::Point const*> stations;
		getStationsInsideSubPolygons(stations);
		dynamic_cast<GMSHAdaptiveMeshDensity*>(_mesh_density_strategy)->addPoints(stations);
	}
}

void GMSHPolygonTree::createGMSHPoints(std::vector<FileIO::GMSHPoint*> & gmsh_pnts) const
{
	const size_t n_pnts_polygon (_node_polygon->getNumberOfPoints());
	for (size_t k(0); k<n_pnts_polygon; k++) {
		const size_t id (_node_polygon->getPointID(k));
		GEOLIB::Point const*const pnt(_node_polygon->getPoint(k));
		gmsh_pnts[id] = new GMSHPoint(*pnt, id, _mesh_density_strategy->getMeshDensityAtPoint(pnt));
	}

	const size_t n_plys(_plys.size());
	for (size_t k(0); k<n_plys; k++) {
		const size_t n_pnts_in_ply(_plys[k]->getNumberOfPoints());
		for (size_t j(0); j<n_pnts_in_ply; j++) {
			if (_node_polygon->isPntInPolygon(*(_plys[k]->getPoint(j)))) {
				const size_t id (_plys[k]->getPointID(j));
				GEOLIB::Point const*const pnt(_plys[k]->getPoint(j));
				gmsh_pnts[id] = new GMSHPoint(*pnt, id, _mesh_density_strategy->getMeshDensityAtPoint(pnt));
			}
		}
	}

	// walk through childs
	for (std::list<SimplePolygonTree*>::const_iterator it (_childs.begin()); it != _childs.end(); it++) {
		dynamic_cast<GMSHPolygonTree*>((*it))->createGMSHPoints(gmsh_pnts);
	}
}

void GMSHPolygonTree::writeLineLoop(size_t &line_offset, size_t &sfc_offset, std::ostream& out) const
{
	const size_t n_pnts (_node_polygon->getNumberOfPoints());
	size_t first_pnt_id(_node_polygon->getPointID(0)), second_pnt_id;
	for (size_t k(1); k<n_pnts; k++) {
		second_pnt_id = _node_polygon->getPointID(k);
		if (265 <= line_offset+k-1 && line_offset+k-1 <= 270) {
			std::cout << "Line(" << line_offset + k-1 << ") = {" << first_pnt_id << "," << second_pnt_id << "};" << std::endl;
		}
		out << "Line(" << line_offset + k-1 << ") = {" << first_pnt_id << "," << second_pnt_id << "};" << std::endl;
		first_pnt_id = second_pnt_id;
	}
	out << "Line Loop(" << line_offset + n_pnts-1 << ") = {";
	for (size_t k(0); k<n_pnts - 2; k++) {
		out << line_offset+k << ",";
	}
	out << line_offset+n_pnts-2 << "};" << std::endl;
	out << "Plane Surface(" << sfc_offset << ") = {" << line_offset+n_pnts-1 << "};" << std::endl;
	line_offset += n_pnts;
	sfc_offset++;
}

void GMSHPolygonTree::writeLineConstraints(size_t &line_offset, size_t sfc_number, std::ostream& out) const
{
	const size_t n_plys (_plys.size());
	for (size_t j(0); j<n_plys; j++) {
		const size_t n_pnts(_plys[j]->getNumberOfPoints());
		size_t first_pnt_id(_plys[j]->getPointID(0)), second_pnt_id;
		for (size_t k(1); k<n_pnts; k++) {
			second_pnt_id = _plys[j]->getPointID(k);
			if (_plys[j]->isSegmentMarked(k-1) && _node_polygon->isPntInPolygon(*(_plys[j]->getPoint(k)))) {
				out << "Line(" << line_offset + k-1 << ") = {" << first_pnt_id << "," << second_pnt_id << "};" << std::endl;
				out << "Line { " << line_offset+k-1 << " } In Surface { " << sfc_number << " };" << std::endl;
			}
			first_pnt_id = second_pnt_id;
		}
		line_offset += n_pnts;
	}
}

void GMSHPolygonTree::writeSubPolygonsAsLineConstraints(size_t &line_offset, size_t sfc_number, std::ostream& out) const
{
	for (std::list<SimplePolygonTree*>::const_iterator it (_childs.begin()); it != _childs.end(); it++) {
		dynamic_cast<GMSHPolygonTree*>((*it))->writeSubPolygonsAsLineConstraints(line_offset, sfc_number, out);
	}

	if (_parent != NULL) {
		const size_t n_pnts(_node_polygon->getNumberOfPoints());
		size_t first_pnt_id(_node_polygon->getPointID(0)), second_pnt_id;
		for (size_t k(1); k<n_pnts; k++) {
			second_pnt_id = _node_polygon->getPointID(k);
			out << "Line(" << line_offset + k-1 << ") = {" << first_pnt_id << "," << second_pnt_id << "};" << std::endl;
			first_pnt_id = second_pnt_id;
			out << "Line { " << line_offset+k-1 << " } In Surface { " << sfc_number << " };" << std::endl;
		}
		line_offset += n_pnts;
	}

}
void GMSHPolygonTree::writeStations(size_t & pnt_id_offset, size_t sfc_number, std::ostream& out) const
{
	const size_t n_stations(_stations.size());
	for (size_t k(0); k<n_stations; k++) {
		out << "Point(" << pnt_id_offset + k << ") = {" << (*(_stations[k]))[0] << "," << (*(_stations[k]))[1] << ", 0.0, ";
		out << _mesh_density_strategy->getMeshDensityAtPoint(_stations[k]) << "};" << std::endl;
		out << "Point { " << pnt_id_offset + k << " } In Surface { " << sfc_number << " }; " << std::endl;
	}
	pnt_id_offset += n_stations;
}

void GMSHPolygonTree::getPointsFromSubPolygons(std::vector<GEOLIB::Point const*>& pnts)
{
	const size_t n_pnts_polygon (_node_polygon->getNumberOfPoints());
	for (size_t k(0); k<n_pnts_polygon; k++) {
		pnts.push_back(_node_polygon->getPoint(k));
	}

	for (std::list<SimplePolygonTree*>::const_iterator it (_childs.begin()); it != _childs.end(); it++) {
		dynamic_cast<GMSHPolygonTree*>((*it))->getPointsFromSubPolygons(pnts);
	}
}

void GMSHPolygonTree::getStationsInsideSubPolygons(std::vector<GEOLIB::Point const*>& stations)
{
	const size_t n_stations(_stations.size());
	for (size_t k(0); k<n_stations; k++) {
		stations.push_back(_stations[k]);
	}

	for (std::list<SimplePolygonTree*>::const_iterator it (_childs.begin()); it != _childs.end(); it++) {
		dynamic_cast<GMSHPolygonTree*>((*it))->getStationsInsideSubPolygons(stations);
	}
}

} // end namespace FileIO
