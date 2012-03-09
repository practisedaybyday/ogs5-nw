/*
 * GMSHInterface.cpp
 *
 *  Created on: Apr 29, 2010
 *      Author: TF
 */

#include <fstream>
#include <list>
#include <vector>

// Base
#include "swap.h"
#include "Configure.h"
#include "BuildInfo.h"

// FileIO
#include "GMSHInterface.h"
#include "GMSHNoMeshDensity.h"
#include "GMSHNoMeshDensity.h"
#include "GMSHFixedMeshDensity.h"
#include "GMSHAdaptiveMeshDensity.h"

// GEOLIB
#include "Point.h"
#include "Polygon.h"
#include "Polyline.h"
#include "QuadTree.h"
#include "SimplePolygonTree.h"

// MSH
#include "msh_elem.h"
#include "msh_mesh.h"

namespace FileIO
{
GMSHInterface::GMSHInterface(GEOLIB::GEOObjects const& geo_objs,
				bool include_stations_as_constraints,
				MeshDensityAlgorithm mesh_density_algorithm,
				double param1, double param2, size_t param3,
				std::vector<std::string>& selected_geometries) :
	_n_pnt_offset (0),
	_n_lines (0),
	_n_plane_sfc(0),
	_geo_objs(geo_objs),
	_selected_geometries(selected_geometries),
	_include_stations_as_constraints(include_stations_as_constraints)
{
	if (mesh_density_algorithm == GMSHInterface::NoMeshDensity) {
		_mesh_density_strategy = new GMSHNoMeshDensity;
	} else {
		if(mesh_density_algorithm == GMSHInterface::FixedMeshDensity) {
			_mesh_density_strategy = new GMSHFixedMeshDensity(param1);
		} else {
			_mesh_density_strategy = new GMSHAdaptiveMeshDensity(param1, param2, param3);
		}
	}
}

void GMSHInterface::setSelectedGeometries (std::vector<std::string>& selected_geometries)
{
	_selected_geometries = selected_geometries;
	_n_pnt_offset = 0;
	_n_lines = 0;
	_n_plane_sfc = 0;
}

int GMSHInterface::write(std::ostream& out)
{
	out << "// GMSH input file created by OpenGeoSys " << OGS_VERSION << " build " << BUILD_TIMESTAMP << std::endl;
	out << std::endl;

	writeGMSHInputFile(out);
//	computeBoundingPolygons(_bp_vec);
//	for (std::vector<GEOLIB::Polygon*>::const_iterator it_bp(_bp_vec.begin()); it_bp != _bp_vec.end(); it_bp++) {
//
//	}
}

void GMSHInterface::writeGMSHInputFile(std::ostream& out)
{
#ifndef NDEBUG
	std::cerr << "[GMSHInterface::writeGMSHInputFile] get data from instance of class GEOObjects ... " << std::flush;
#endif
	// get data from geo
	std::vector<GEOLIB::Point*> all_points;
	std::vector<GEOLIB::Polyline*> all_polylines;
	std::vector<GEOLIB::Point*> all_stations;
	fetchGeometries (all_points, all_polylines, all_stations);
#ifndef NDEBUG
	std::cerr << "ok" << std::endl;
#endif
	_mesh_density_strategy->init(all_points);
	if (dynamic_cast<GMSHAdaptiveMeshDensity*> (_mesh_density_strategy)
					&& _include_stations_as_constraints) {
		dynamic_cast<GMSHAdaptiveMeshDensity*> (_mesh_density_strategy)->addPoints(all_stations);
	}

	// search bounding polygon
	size_t bp_idx(0); // bounding polygon index
	GEOLIB::Polygon* bounding_polygon(getBoundingPolygon(all_polylines, bp_idx));
	if (!bounding_polygon)
		return;

	// *** GMSH - write all non-station points
	size_t pnts_written = writeGMSHPointsInsideBoundingPolygon(bounding_polygon, all_points, out, false);

	// write the bounding polygon
	writeBoundingPolygon(bounding_polygon, _n_pnt_offset-pnts_written, out);

	// write all other polylines as constraints
	writePolylinesAsConstraints(bounding_polygon, bp_idx, all_polylines, out);

	// write stations as constraints
	out << "// Stations" << std::endl;
	writeGMSHPoints(all_stations, out, true);

	// write Steiner points
	if (dynamic_cast<GMSHAdaptiveMeshDensity*> (_mesh_density_strategy)) {
		out << "// Steiner points" << std::endl;
		dynamic_cast<GMSHAdaptiveMeshDensity*> (_mesh_density_strategy)->writeSteinerPoints(bounding_polygon, _n_pnt_offset, _n_plane_sfc, out);
	}
}

void GMSHInterface::writeGMSHPoints(const std::vector<GEOLIB::Point*> &pnt_vec, std::ostream& out, bool station)
{
	// write points
	const size_t n (pnt_vec.size());
	if (!station || dynamic_cast<GMSHAdaptiveMeshDensity*>(_mesh_density_strategy) == NULL) {
		for (size_t k(0); k < n; k++) {
			out << "Point(" << _n_pnt_offset + k << ") = {" << (*(pnt_vec[k]))[0] << ","
				 << (*(pnt_vec[k]))[1] << "," << (*(pnt_vec[k]))[2];
			_mesh_density_strategy->getMeshDensityAtPoint(pnt_vec[k], out);
			out << "};" << std::endl;
		}
	} else {
		for (size_t k(0); k < n; k++) {
			out << "Point(" << _n_pnt_offset + k << ") = {" << (*(pnt_vec[k]))[0] << ","
				 << (*(pnt_vec[k]))[1] << "," << (*(pnt_vec[k]))[2];
			dynamic_cast<GMSHAdaptiveMeshDensity*>(_mesh_density_strategy)->getMeshDensityAtStation(pnt_vec[k], out);
			out << "};" << std::endl;
		}
	}
	_n_pnt_offset += n;
}

size_t GMSHInterface::writeGMSHPointsInsideBoundingPolygon(GEOLIB::Polygon const*const bounding_polygon, const std::vector<GEOLIB::Point*> &pnt_vec, std::ostream& out, bool station)
{
	// write points
	const size_t n (pnt_vec.size());
	size_t cnt (0);
	if (!station || dynamic_cast<GMSHAdaptiveMeshDensity*>(_mesh_density_strategy) == NULL) {
		for (size_t k(0); k < n; k++) {
			if (bounding_polygon->isPntInPolygon(*(pnt_vec[k]))) {
				out << "Point(" << _n_pnt_offset + cnt << ") = {" << (*(pnt_vec[k]))[0] << ","
					 << (*(pnt_vec[k]))[1] << "," << (*(pnt_vec[k]))[2];
				_mesh_density_strategy->getMeshDensityAtPoint(pnt_vec[k], out);
				out << "};" << std::endl;
				cnt++;
			}
		}
	} else {
		for (size_t k(0); k < n; k++) {
			if (bounding_polygon->isPntInPolygon(*(pnt_vec[k]))) {
				out << "Point(" << _n_pnt_offset + cnt << ") = {" << (*(pnt_vec[k]))[0] << ","
					 << (*(pnt_vec[k]))[1] << "," << (*(pnt_vec[k]))[2];
				dynamic_cast<GMSHAdaptiveMeshDensity*>(_mesh_density_strategy)->getMeshDensityAtStation(pnt_vec[k], out);
				out << "};" << std::endl;
				cnt++;
			}
		}
	}
	_n_pnt_offset += cnt;
	return cnt;
}


void GMSHInterface::writeGMSHPolyline(std::vector<GEOLIB::Polyline*> const& ply_vec, size_t ply_id,
				size_t offset, std::ostream& out)
{
	size_t local_offset (this->_n_pnt_offset - offset);
	GEOLIB::Polyline const*const ply (ply_vec[ply_id]);
	size_t s (ply->getNumberOfPoints());

	_ply_id_gmsh_line_mapping.push_back(PolylineGMSHMapping(ply_id, _n_lines, _n_lines+s-2));

	// write line segments (= Line) of the polyline
	for (size_t j(0); j < s - 1; j++) {
		out << "Line(" << _n_lines + j << ") = {" << local_offset + ply->getPointID(j) << ","
						<< local_offset + ply->getPointID(j + 1) << "};" << std::endl;
	}
	// write the line segments contained in the polyline (=Line Loop)
	out << "Line Loop (" << _n_lines + s - 1 << ") = {";
	for (size_t j(0); j < s - 2; j++)
		out << _n_lines + j << ",";
	out << _n_lines + s - 2 << "};" << std::endl;
	_n_lines += s;
}

void GMSHInterface::writeGMSHPolylines(const std::vector<GEOLIB::Polyline*>& ply_vec, std::ostream &out)
{
	size_t n (ply_vec.size());
	for (size_t k(0); k < n; k++) {
		// write k-th polyline
		writeGMSHPolyline (ply_vec, k, 0, out);
	}
	out << std::endl;
}

size_t GMSHInterface::writeGMSHPolygon(const std::vector<GEOLIB::Polyline*>& ply_vec, size_t id, size_t offset, std::ostream & out)
{
	writeGMSHPolyline (ply_vec, id, offset, out);
	return _n_lines - 1;
}

//bool GMSHInterface::writeGMSHInputFile(const std::string &proj_name, const GEOLIB::GEOObjects& geo,
//				bool useStationsAsContraints, bool useSteinerPoints, std::ostream & out)
//{
//	std::cerr << "GMSHInterface::writeGMSHInputFile " << std::endl;
//	std::cerr << "get data from geo ... " << std::flush;
//	// get data from geo
//	const std::vector<GEOLIB::Point*>* pnts(geo.getPointVec(proj_name));
//	const std::vector<GEOLIB::Polyline*>* plys(geo.getPolylineVec(proj_name));
//	std::cerr << "ok" << std::endl;
//
//	std::vector<GEOLIB::Point*> station_points;
//	if (useStationsAsContraints) {
//		getStationPoints(geo, station_points);
//		if (dynamic_cast<GMSHAdaptiveMeshDensity*> (_mesh_density_strategy)) {
//			dynamic_cast<GMSHAdaptiveMeshDensity*> (_mesh_density_strategy)->addPoints(station_pnts);
//		}
//	}
//
//	// write points
//	writeGMSHPoints(*pnts, out, false);
//
//	// write Polylines
//	std::map<size_t, size_t> geo2gmsh_polygon_id_map;
//	std::map<size_t, size_t> geo2gmsh_surface_id_map;
//
//	for (size_t i = 0; i < plys->size(); i++) {
//		if ((*plys)[i]->isClosed()) {
//			GEOLIB::Polygon const& polygon(*((*plys)[i]));
//			// omit creating polygons twice
//			size_t k(0);
//			for (k = 0; k < i; k++) {
//				GEOLIB::Polygon const& test_polygon(*((*plys)[k]));
//				if (polygon == test_polygon) {
//					k = i + 1;
//				}
//			}
//			// k==i means polygon is not yet written
//			if (k == i) {
//				size_t polygon_id = writeGMSHPolygon(*plys, i, pnts->size(), out);
//				geo2gmsh_polygon_id_map[i] = polygon_id;
//			}
//		} else writeGMSHPolyline(*plys, i, pnts->size(), out);
//	}
//
//	for (size_t i = 0; i < plys->size(); i++)
//		if ((*plys)[i]->isClosed()) {
//			std::list<size_t> polygon_list(findHolesInsidePolygon(plys, i));
//
//			// search for holes that are twice
//			std::list<size_t>::iterator it_i(polygon_list.begin());
//			while (it_i != polygon_list.end()) {
//				std::list<size_t>::iterator it_j(it_i);
//				it_j++;
//				while (it_j != polygon_list.end()) {
//					bool erased(false);
//					GEOLIB::Polygon const& polygon_i(*((*plys)[*it_i]));
//					GEOLIB::Polygon const& polygon_j(*((*plys)[*it_j]));
//					if (polygon_i == polygon_j) {
//						it_j = polygon_list.erase(it_j);
//						erased = true;
//					}
//					if (!erased) it_j++;
//				}
//				it_i++;
//			}
//			// end search for holes that are twice
//
//			// is polygon a hole?
//			size_t j(0);
//			if (!isPolygonInOtherPolygon(plys, i, j)) {
//				std::list<size_t> polygon_list_gmsh;
//				// apply polygon id mapping (ogs geometry -> gmsh)
//				for (std::list<size_t>::iterator it(polygon_list.begin()); it != polygon_list.end(); it++) {
//					polygon_list_gmsh.push_back((geo2gmsh_polygon_id_map.find(*it))->second);
//				}
//
//				bool write_holes(false);
//				this->writePlaneSurface(polygon_list_gmsh, write_holes);
//				geo2gmsh_surface_id_map[i] = _n_plane_sfc - 1;
//				// write holes as constraints
//				if (!write_holes && polygon_list.size() > 1) {
//					std::list<size_t>::const_iterator it(polygon_list.begin());
//					for (it++; it != polygon_list.end(); it++) {
//						_out << "// adding constraints" << std::endl;
//
//						// search polygon in _ply_id_gmsh_line_mapping
//						size_t k(0);
//						const size_t
//										ply_id_gmsh_line_mapping_size(
//														_ply_id_gmsh_line_mapping.size());
//						for (; k < ply_id_gmsh_line_mapping_size; k++) {
//							if (_ply_id_gmsh_line_mapping[k]._ply_id == *it) {
//								for (size_t l(_ply_id_gmsh_line_mapping[k]._gmsh_line_start_id); l
//												< _ply_id_gmsh_line_mapping[k]._gmsh_line_end_id; l++) {
//									_out << "Line {" << l << "} In Surface {" << _n_plane_sfc - 1
//													<< "};" << std::endl;
//								}
//								break;
//							}
//						}
//					}
//				}
//			}
//		}
//
//	if (useStationsAsContraints) addPointsAsConstraints(station_points, *plys,
//					geo2gmsh_surface_id_map, out);
//
//	if (useSteinerPoints) {
//		std::vector<GEOLIB::Point*> steiner_points;
//		getSteinerPoints(quad_tree, steiner_points);
//		_out << "// Steiner points" << std::endl;
//		addPointsAsConstraints(steiner_points, *plys, geo2gmsh_surface_id_map, quad_tree, 0.3);
//		for (size_t i = 0; i < steiner_points.size(); i++)
//			delete steiner_points[i];
//	}
//
//	delete quad_tree;
//
//	std::cerr << "ok" << std::endl;
//
//	return true;
//}

std::list<size_t> GMSHInterface::findHolesInsidePolygon(const std::vector<GEOLIB::Polyline*>* plys,
				size_t i) const
{
	GEOLIB::Polygon polygon(*((*plys)[i]));
	std::list<size_t> polygon_list;
	polygon_list.push_back(i);
	for (size_t j = 0; j < plys->size(); j++) {
		// check if polygons are located completely inside the given polygon
		if ((i != j) && ((*plys)[j]->isClosed()))
		{
			GEOLIB::Polyline* line ((*plys)[j]);
			bool isInside(true);
			for (size_t k = 0; k < line->getNumberOfPoints(); k++)
				if (!polygon.isPntInPolygon(*(line->getPoint(k))))
				{
					isInside = false;
					break;
				}
			if (isInside)
				polygon_list.push_back(j);
		}
	}
	return polygon_list;
}

bool GMSHInterface::isPolygonInOtherPolygon(const std::vector<GEOLIB::Polyline*>* plys, size_t i, size_t &j) const
{
	if (! ((*plys)[i])->isClosed())
		return false;

	GEOLIB::Polygon const& ply_i(*((*plys)[i])); // polygon to test with
	const size_t n_pnts_ply_i(ply_i.getNumberOfPoints());
	const size_t n_plys(plys->size());
	bool is_not_hole(true);
	// check other polygons
	for (j = 0; j < n_plys && is_not_hole; j++) {
		if ((i != j) && (((*plys)[j])->isClosed())) {
			GEOLIB::Polygon const& ply_j(*((*plys)[j]));
			bool is_hole_in_ply_j(true);
			for (size_t k(0); k < n_pnts_ply_i && is_hole_in_ply_j; k++) {
				if (!ply_j.isPntInPolygon(*(ply_i.getPoint(k)))) {
					is_hole_in_ply_j = false;
				}
			}
			if (is_hole_in_ply_j) {
				is_not_hole = false;
			}
		}
	}
	if (i==j) return false;
	return !is_not_hole;
}

void GMSHInterface::writePlaneSurface (std::list<size_t> const& polygon_list, bool respect_holes)
{
	_out << "Plane Surface (" << _n_plane_sfc << ") = {" << std::flush;
	std::list<size_t>::const_iterator it (polygon_list.begin());
	_out << *it << std::flush;
	if (respect_holes) {
		for (++it; it != polygon_list.end(); ++it)
			_out << ", " << *it << std::flush;
	}
	_out << "};" << std::endl;
	_n_plane_sfc++;
}



//void GMSHInterface::writeAllDataToGMSHInputFileAdaptive (
//        GEOLIB::GEOObjects& geo,
//        std::vector<std::string> const &selected_geometries,
//        size_t number_of_point_per_quadtree_node,
//        double mesh_density_scaling,
//        double mesh_density_scaling_station_pnts)
//{
//	std::string file_string = this->writeAllDataToGMSHInputAdaptive(geo, selected_geometries, number_of_point_per_quadtree_node, mesh_density_scaling, mesh_density_scaling_station_pnts);
//}

//void GMSHInterface::writeAllDataToGMSHInputFileFixedMeshDensity (std::ostream &out,
//        double mesh_density)
//{
//#ifndef NDEBUG
//	std::cout << "[GMSHInterface::writeGMSHInputFile] non adaptive" << std::endl;
//#endif
//
//	std::vector<GEOLIB::Point*> all_points;
//	std::vector<GEOLIB::Polyline*> all_polylines;
//	std::vector<GEOLIB::Point*> all_stations;
//	fetchGeometries (all_points, all_polylines, all_stations);
//
//	// search bounding polygon
//	size_t bp_idx (0); // bounding polygon index
//	GEOLIB::Polygon* bounding_polygon (getBoundingPolygon(all_polylines, bp_idx));
//	if (!bounding_polygon)
//		return;
//
//	// *** GMSH - write all non-station points
//	writeGMSHPoints(all_points, out, false);
//
//#ifndef NDEBUG
//	std::cout << "\twrite bounding polygon ... " << std::flush;
//#endif
//	// write bounding polygon
//	writeBoundingPolygon (bounding_polygon, _n_pnt_offset - all_points.size(), out);
//
//	// write all other polylines as constraints
//	writePolylinesAsConstraints(bounding_polygon, bp_idx, all_polylines, mesh_density, NULL);
//
//	// write stations as constraints
//	out << "// Stations" << std::endl;
//	writeGMSHPoints(all_stations, out, true);
//#ifndef NDEBUG
//	std::cout << "ok" << std::endl;
//#endif
//}

void GMSHInterface::fetchGeometries (std::vector<GEOLIB::Point*>& all_points,
                                     std::vector<GEOLIB::Polyline*>& all_polylines,
                                     std::vector<GEOLIB::Point*>& all_stations) const
{
	// get names of all available data sources except stations
	std::vector<std::string> geo_names;
	// get station names
	std::vector<std::string> geo_station_names;

	for (std::vector<std::string>::const_iterator it (_selected_geometries.begin());
	     it != _selected_geometries.end(); ++it)
	{
		if ((_geo_objs.getPointVecObj (*it))->getType() == GEOLIB::PointVec::POINT)
			geo_names.push_back (*it);
		else if ((_geo_objs.getPointVecObj (*it))->getType() == GEOLIB::PointVec::STATION)
			geo_station_names.push_back (*it);
	}

	size_t pnt_offset (0);
	// fetch points and polylines and add them to the vectors, add points to the QuadTree
	for (std::vector<std::string>::const_iterator it (geo_names.begin());
	     it != geo_names.end(); ++it)
	{
		// get data from geo
#ifndef NDEBUG
		std::cout << "\tfetch geometrical data for " << *it << " " << std::flush;
#endif
		const std::vector<GEOLIB::Point*>* pnts (_geo_objs.getPointVec (*it));
		const std::vector<GEOLIB::Polyline*>* plys (_geo_objs.getPolylineVec (*it));
#ifndef NDEBUG
		std::cout << "ok" << std::endl;
#endif

		if (pnts) {
			// insert points into vector all_points
			all_points.insert (all_points.end(), pnts->begin(), pnts->end());

			if (plys) {
				for (size_t k(0); k < plys->size(); k++) {
					size_t pos (all_polylines.size());
					// insert new polyline
					all_polylines.push_back (new GEOLIB::Polyline (all_points));
					// copy points
					for (size_t j(0); j < (*plys)[k]->getNumberOfPoints(); j++) {
						// set points of polyline
						(all_polylines[pos])->addPoint(pnt_offset + ((*plys)[k])-> getPointID(j));
					}
				}
			}
			pnt_offset += pnts->size();
		}
	}

	for (std::vector<std::string>::const_iterator it (geo_station_names.begin());
	     it != geo_station_names.end(); ++it)
	{
		// get data from geo
#ifndef NDEBUG
		std::cout << "\tfetch station data for " << *it << " " << std::flush;
#endif
		const std::vector<GEOLIB::Point*>* pnts (_geo_objs.getPointVec (*it));
#ifndef NDEBUG
		std::cout << "ok" << std::endl;
#endif
		// insert points into vector all_stations
		all_stations.insert (all_stations.end(), pnts->begin(), pnts->end());
	}
}

GEOLIB::Polygon* GMSHInterface::getBoundingPolygon(
				std::vector<GEOLIB::Polyline*> const & all_polylines, size_t& bp_idx) const
{
	GEOLIB::Polygon* bounding_polygon(NULL);
	const size_t n_polylines(all_polylines.size());
	for (size_t k(0); k < n_polylines; k++) {
		if (all_polylines[k]->isClosed()) // == Polygon
		{
			if (bounding_polygon) // we have already a bounding polygon
			{
				if (!bounding_polygon->isPolylineInPolygon(*(all_polylines[k]))) {
					GEOLIB::Polygon* tmp_polygon(new GEOLIB::Polygon(*(all_polylines[k])));
					if (tmp_polygon->isPolylineInPolygon(*bounding_polygon)) {
						// found new bounding polygon
						delete bounding_polygon;
						bounding_polygon = tmp_polygon;
						bp_idx = k;
					} else std::cerr
									<< "INFO: there is no inclusion relation between the polygons "
									<< k << " and " << bp_idx << std::endl;
				}
			} else {
				bounding_polygon = new GEOLIB::Polygon(*(all_polylines[k]));
				bp_idx = k;
			}
		}
	}

	if (!bounding_polygon) {
		std::cerr
						<< "WARNING: GMSHInterface::writeAllDataToGMSHInputFile: did not found bounding polygon - abort writing"
						<< std::endl;
		return NULL;
	}

	return bounding_polygon;
}

void GMSHInterface::writeBoundingPolygon(GEOLIB::Polygon const* const bounding_polygon,
				size_t pnt_offset, std::ostream & out)
{
	std::cout << "write bounding polygon ... " << std::flush;
	// write bounding polygon
	size_t s (bounding_polygon->getNumberOfPoints());
	// write line segments (= Line) of the polyline
	for (size_t j(0); j < s - 1; j++)
		out << "Line(" << _n_lines + j << ") = {" <<  pnt_offset + bounding_polygon->getPointID(j) << ","
		     << pnt_offset + bounding_polygon->getPointID(j + 1) << "};" << std::endl;
	// write the line segments contained in the polyline (=Line Loop)
	out << "Line Loop (" << _n_lines + s - 1 << ") = {";
	for (size_t j(0); j < s - 2; j++)
		out << _n_lines + j << ",";
	out << _n_lines + s - 2 << "};" << std::endl;
	_n_lines += s;
	// write plane surface
	out << "Plane Surface (" << _n_plane_sfc << ") = {" << _n_lines - 1 << "};" << std::endl;
	_n_plane_sfc++;
	std::cout << "ok" << std::endl;
}

//void GMSHInterface::addPointsAsConstraints(const std::vector<GEOLIB::Point*> &points,
//				const std::vector<GEOLIB::Polyline*> &polylines, std::ostream &out)
//{
//	std::vector<GEOLIB::Polygon*> polygons;
//	for (size_t j = 0; j < polylines.size(); j++) {
//		if (polylines[j]->isClosed()) {
//			GEOLIB::Polygon* pgn = new GEOLIB::Polygon(*polylines[j]);
//			polygons.push_back(pgn);
//		}
//	}
//
//	size_t nPoints = points.size();
//	for (size_t i = 0; i < nPoints; i++) {
//		std::list<size_t> surrounding_polygons;
//		for (size_t j = 0; j < polygons.size(); j++)
//			if (polygons[j]->isPntInPolygon(*(points[i])))
//				surrounding_polygons.push_back(j);
//
//		if (!surrounding_polygons.empty()) {
//			for (std::list<size_t>::iterator it = surrounding_polygons.begin(); it
//							!= surrounding_polygons.end(); ++it) {
//				for (std::list<size_t>::iterator jt = surrounding_polygons.begin(); jt
//								!= surrounding_polygons.end();) {
//					if (it != jt) {
//						if (polygons[*it]->isPolylineInPolygon(*(polygons[*jt])))
//							it = surrounding_polygons.erase(it);
//						else ++it;
//					} else ++jt;
//				}
//			}
//			_n_pnt_offset++;
//			out << "Point(" << _n_pnt_offset << ") = {" << (*points[i])[0] << ","
//				<< (*points[i])[1] << "," << (*points[i])[2]
//				<< _mesh_density_strategy->getMeshDensityAtPoint(points[i], out) << "};"
//							<< std::endl;
//			out << "Point {" << _n_pnt_offset << "} In Surface {"
//							<< geo2gmsh_surface_id_map[*(surrounding_polygons.begin())] << "};"
//							<< std::endl;
//		}
//	}
//}

void GMSHInterface::getStationPoints(const GEOLIB::GEOObjects &geo, std::vector<GEOLIB::Point*>& station_points) const
{
	std::vector<std::string> stn_names;
	geo.getStationVectorNames(stn_names);

	for (std::vector<std::string>::const_iterator it (stn_names.begin()); it != stn_names.end();
	     ++it)
	{
		const std::vector<GEOLIB::Point*>* pnts (geo.getPointVec (*it));
		station_points.insert (station_points.end(), pnts->begin(), pnts->end());
	}
}

//void GMSHInterface::getSteinerPoints(
//        GEOLIB::QuadTree<GEOLIB::Point>* quad_tree,
//        std::vector<GEOLIB::Point*>& steiner_points) const
//{
//	std::list<GEOLIB::QuadTree<GEOLIB::Point>*> leaf_list;
//	quad_tree->getLeafs (leaf_list);
//	for (std::list<GEOLIB::QuadTree<GEOLIB::Point>*>::const_iterator it (leaf_list.begin());
//	     it != leaf_list.end(); it++)
//		if ((*it)->getPoints().empty())
//		{
//			// compute point from square
//			GEOLIB::Point ll, rr;
//			(*it)->getSquarePoints (ll, rr);
//			GEOLIB::Point* mid_point =
//			        new GEOLIB::Point(0.5 * (rr[0] + ll[0]),
//			                          0.5 * (rr[1] + ll[1]),
//			                          0.5 * (rr[2] + ll[2]));
//			steiner_points.push_back(mid_point);
//		}
//}

void GMSHInterface::writePolylinesAsConstraints(GEOLIB::Polygon * bounding_polygon, size_t bp_idx,
				std::vector<GEOLIB::Polyline*> const& all_polylines, std::ostream & out)
{
	const size_t n_polylines(all_polylines.size());
	for (size_t k(0); k < n_polylines; k++) {
		if (k != bp_idx) {
			bool begin_line_pnt_inside_polygon(true);
			bool end_line_pnt_inside_polygon(true);

			size_t s(all_polylines[k]->getNumberOfPoints());

			// write line segments (= Line) of the polyline
			for (size_t j(0); j < s - 1; j++) {
				// check if line segment is contained in bounding polygon
				bool line_seg_is_already_used(GEOLIB::containsEdge(
								*(dynamic_cast<GEOLIB::Polyline*> (bounding_polygon)),
								(all_polylines[k])-> getPointID(j),
								(all_polylines[k])-> getPointID(j + 1)));
				// check if line segment is contained in a previous polyline
				for (size_t i(0); i < k && !line_seg_is_already_used; i++)
					line_seg_is_already_used = GEOLIB::containsEdge(*(all_polylines[i]),
									(all_polylines[k])-> getPointID(j),
									(all_polylines[k])->getPointID(j + 1));

				if (!line_seg_is_already_used) {
					// check if first point of polyline is inside bounding polygon
					if (j == 0) begin_line_pnt_inside_polygon = bounding_polygon->isPntInPolygon(
									*(all_polylines[k])->getPoint(j));
					// check if end point of the line is inside bounding polygon
					end_line_pnt_inside_polygon = bounding_polygon->isPntInPolygon(
									*(all_polylines[k])->getPoint(j + 1));

					if (begin_line_pnt_inside_polygon && end_line_pnt_inside_polygon) {
						out << "Line(" << _n_lines + j << ") = {"
										<< (all_polylines[k])->getPointID(j) << ","
										<< (all_polylines[k])->getPointID(j + 1) << "};"
										<< std::endl;
						// write line as constraint
						out << "Line {" << _n_lines + j << "} In Surface {" << _n_plane_sfc - 1
										<< "};" << std::endl;
					} else {
						if (begin_line_pnt_inside_polygon && !end_line_pnt_inside_polygon) {
							// create new point
							GEOLIB::Point * s(bounding_polygon-> getIntersectionPointPolygonLine(
											*(all_polylines[k])-> getPoint(j),
											*(all_polylines[k])-> getPoint(j + 1)));
							if (s != NULL) {
								out << "Point(" << _n_pnt_offset << ") = {" << (*s)[0] << ","
												<< (*s)[1] << "," << (*s)[2];
								_mesh_density_strategy->getMeshDensityAtPoint(all_polylines[k]->getPoint(j), out);
								out << "}; // new end point of polyline " << k << std::endl;
								// write line
								out << "Line(" << _n_lines + j << ") = {"
												<< (all_polylines[k])-> getPointID(j) << ","
												<< _n_pnt_offset << "};" << std::endl;
								// write line as constraint
								out << "Line {" << _n_lines + j << "} In Surface {"
												<< _n_plane_sfc - 1 << "};" << std::endl;
								_n_pnt_offset++;
								delete s;
							}
						}
						if (!begin_line_pnt_inside_polygon && end_line_pnt_inside_polygon) {
							// create new point
							GEOLIB::Point * s(bounding_polygon-> getIntersectionPointPolygonLine(
											*(all_polylines[k])-> getPoint(j),
											*(all_polylines[k])-> getPoint(j + 1)));
							if (s != NULL) {
								out << "Point(" << _n_pnt_offset << ") = {" << (*s)[0] << "," << (*s)[1] << "," << (*s)[2];
								_mesh_density_strategy->getMeshDensityAtPoint(all_polylines[k]->getPoint(j), out);
								out << "}; // new end point of polyline " << k << std::endl;
								// write line
								out << "Line(" << _n_lines + j << ") = {" << _n_pnt_offset << ","
												<< (all_polylines[k])->getPointID(j + 1) << "};"
												<< std::endl;
								// write line as constraint
								out << "Line {" << _n_lines + j << "} In Surface {"
												<< _n_plane_sfc - 1 << "};" << std::endl;
								_n_pnt_offset++;
								delete s;
							}
						}
					}
					begin_line_pnt_inside_polygon = end_line_pnt_inside_polygon;
				}
			}
			// update line counter
			_n_lines += s;
		}
	}
}

bool GMSHInterface::isGMSHMeshFile (const std::string& fname)
{
	std::ifstream input (fname.c_str());

	if (!input)
	{
		std::cerr << "GMSHInterface::isGMSHMeshFile could not open file " << fname <<
		std::endl;
		return false;
	}

	std::string header_first_line;
	input >> header_first_line;
	if (header_first_line.find ("$MeshFormat") != std::string::npos)
	{
		// read version
		std::string version;
		getline (input, version);
		getline (input, version);
		std::cerr << "found GMSH mesh file version: " << version << std::endl;
		input.close ();
		return true;
	}

	return false;
}

void GMSHInterface::readGMSHMesh(std::string const& fname,
                                 MeshLib::CFEMesh* mesh)
{
	std::string line;
	std::ifstream in (fname.c_str(), std::ios::in);
	getline(in, line); // Node keyword

	if (line.find("$MeshFormat") != std::string::npos)
	{
		getline(in, line); // version-number file-type data-size
		getline(in, line); //$EndMeshFormat
		getline(in, line); //$Nodes Keywords

		size_t n_nodes(0);
		size_t n_elements(0);
		while (line.find("$EndElements") == std::string::npos)
		{
			// Node data
			long id;
			double x, y, z;
			in >> n_nodes >> std::ws;
			for (size_t i = 0; i < n_nodes; i++)
			{
				in >> id >> x >> y >> z >> std::ws;
				mesh->nod_vector.push_back(new MeshLib::CNode(id, x, y, z));
			}
			getline(in, line); // End Node keyword $EndNodes

			// Element data
			getline(in, line); // Element keyword $Elements
			in >> n_elements >> std::ws; // number-of-elements
			for (size_t i = 0; i < n_elements; i++)
			{
				MeshLib::CElem* elem (new MeshLib::CElem(i));
				elem->Read(in, 7);
				if (elem->GetElementType() != MshElemType::INVALID)
					mesh->ele_vector.push_back(elem);
			}
			getline(in, line); // END keyword

			// correct indices TF
			const size_t n_elements(mesh->ele_vector.size());
			for (size_t k(0); k < n_elements; k++)
				mesh->ele_vector[k]->SetIndex(k);

			// ordering nodes and closing gaps TK
			std::vector<size_t> gmsh_id;
			size_t counter(0);
			for (size_t i = 0; i < mesh->nod_vector.size(); i++)
			{
				const size_t diff = mesh->nod_vector[i]->GetIndex() - counter;
				if (diff == 0)
				{
					gmsh_id.push_back(i);
					counter++;
				}
				else
				{
					for (size_t j = 0; j < diff; j++)
					{
						gmsh_id.push_back(i);
						counter++;
					}
					i--;
				}
			}

			for (size_t i = 0; i < mesh->ele_vector.size(); i++)
				for (long j = 0; j < mesh->ele_vector[i]->GetVertexNumber(); j++)
					mesh->ele_vector[i]->getNodeIndices()[j] =
					        gmsh_id[mesh->ele_vector[i]->GetNodeIndex(j) + 1];

			for (size_t i = 0; i < mesh->nod_vector.size(); i++)
				mesh->nod_vector[i]->SetIndex(i);
			// END OF: ordering nodes and closing gaps TK
		} /*End while*/
	}
	in.close();
}
} // end namespace FileIO
