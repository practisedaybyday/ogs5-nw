/*
 * GMSHInterface.cpp
 *
 *  Created on: Apr 29, 2010
 *      Author: TF
 */

#include <fstream>
#include <vector>

// Base
#include "swap.h"
#include "Configure.h"
#include "BuildInfo.h"

// FileIO
#include "GMSHInterface.h"
#include "GMSHNoMeshDensity.h"
#include "GMSHFixedMeshDensity.h"
#include "GMSHAdaptiveMeshDensity.h"

// GEOLIB
#include "Point.h"
#include "Polygon.h"
#include "Polyline.h"
#include "QuadTree.h"
#include "PolylineWithSegmentMarker.h"

// MSH
#include "msh_elem.h"
#include "msh_mesh.h"

namespace FileIO {
GMSHInterface::GMSHInterface(GEOLIB::GEOObjects & geo_objs, bool include_stations_as_constraints,
				GMSH::MeshDensityAlgorithm mesh_density_algorithm, double param1, double param2,
				size_t param3, std::vector<std::string>& selected_geometries) :
	_n_lines(0), _n_plane_sfc(0), _geo_objs(geo_objs), _selected_geometries(selected_geometries),
	_include_stations_as_constraints(include_stations_as_constraints)
{
	switch (mesh_density_algorithm) {
	case GMSH::NoMeshDensity:
		_mesh_density_strategy = new GMSHNoMeshDensity;
		break;
	case GMSH::FixedMeshDensity:
		_mesh_density_strategy = new GMSHFixedMeshDensity(param1);
		break;
	case GMSH::AdaptiveMeshDensity:
		_mesh_density_strategy = new GMSHAdaptiveMeshDensity(param1, param2, param3);
		break;
	}
}

bool GMSHInterface::isGMSHMeshFile(const std::string& fname)
{
	std::ifstream input(fname.c_str());

	if (!input) {
		std::cerr << "GMSHInterface::isGMSHMeshFile could not open file " << fname << std::endl;
		return false;
	}

	std::string header_first_line;
	input >> header_first_line;
	if (header_first_line.find("$MeshFormat") != std::string::npos) {
		// read version
		std::string version;
		getline(input, version);
		getline(input, version);
		std::cerr << "found GMSH mesh file version: " << version << std::endl;
		input.close();
		return true;
	}

	return false;
}

void GMSHInterface::readGMSHMesh(std::string const& fname, MeshLib::CFEMesh* mesh)
{
	std::string line;
	std::ifstream in(fname.c_str(), std::ios::in);
	getline(in, line); // Node keyword

	if (line.find("$MeshFormat") != std::string::npos) {
		getline(in, line); // version-number file-type data-size
		getline(in, line); //$EndMeshFormat
		getline(in, line); //$Nodes Keywords

		size_t n_nodes(0);
		size_t n_elements(0);
		while (line.find("$EndElements") == std::string::npos) {
			// Node data
			long id;
			double x, y, z;
			in >> n_nodes >> std::ws;
			for (size_t i = 0; i < n_nodes; i++) {
				in >> id >> x >> y >> z >> std::ws;
				mesh->nod_vector.push_back(new MeshLib::CNode(id, x, y, z));
			}
			getline(in, line); // End Node keyword $EndNodes

			// Element data
			getline(in, line); // Element keyword $Elements
			in >> n_elements >> std::ws; // number-of-elements
			for (size_t i = 0; i < n_elements; i++) {
				MeshLib::CElem* elem(new MeshLib::CElem(i));
				elem->Read(in, 7);
				if (elem->GetElementType() != MshElemType::INVALID) mesh->ele_vector.push_back(elem);
			}
			getline(in, line); // END keyword

			// correct indices TF
			const size_t n_elements(mesh->ele_vector.size());
			for (size_t k(0); k < n_elements; k++)
				mesh->ele_vector[k]->SetIndex(k);

			// ordering nodes and closing gaps TK
			std::vector<size_t> gmsh_id;
			size_t counter(0);
			for (size_t i = 0; i < mesh->nod_vector.size(); i++) {
				const size_t diff = mesh->nod_vector[i]->GetIndex() - counter;
				if (diff == 0) {
					gmsh_id.push_back(i);
					counter++;
				} else {
					for (size_t j = 0; j < diff; j++) {
						gmsh_id.push_back(i);
						counter++;
					}
					i--;
				}
			}

			for (size_t i = 0; i < mesh->ele_vector.size(); i++)
				for (long j = 0; j < mesh->ele_vector[i]->GetVertexNumber(); j++)
					mesh->ele_vector[i]->getNodeIndices()[j]
									= gmsh_id[mesh->ele_vector[i]->GetNodeIndex(j) + 1];

			for (size_t i = 0; i < mesh->nod_vector.size(); i++)
				mesh->nod_vector[i]->SetIndex(i);
			// END OF: ordering nodes and closing gaps TK
		} /*End while*/
	}
	in.close();
}

int GMSHInterface::write(std::ostream& out)
{
	out << "// GMSH input file created by OpenGeoSys " << OGS_VERSION << " build "
					<< BUILD_TIMESTAMP << std::endl;
	out << std::endl;

	writeGMSHInputFile(out);
	return 1;
}

void GMSHInterface::writeGMSHInputFile(std::ostream& out)
{
#ifndef NDEBUG
	std::cerr << "[GMSHInterface::writeGMSHInputFile] get data from GEOObjects ... " << std::flush;
#endif
	// *** get and merge data from _geo_objs
	_gmsh_geo_name = "GMSHGeometry";
	_gmsh_station_name = "GMSHStations";
	_geo_objs.mergeGeometries(_selected_geometries, _gmsh_geo_name);
	_geo_objs.mergeStations(_selected_geometries, _gmsh_station_name);
	std::vector<GEOLIB::Point*> const* merged_pnts(_geo_objs.getPointVec(_gmsh_geo_name));
	if (! merged_pnts) {
		std::cerr << "[GMSHInterface::writeGMSHInputFile] did not found any points" << std::endl;
		return;
	}
	std::vector<GEOLIB::Polyline*> const* merged_plys(_geo_objs.getPolylineVec(_gmsh_geo_name));
	std::vector<GEOLIB::Point*> const* merged_stations(_geo_objs.getStationVec(_gmsh_station_name));
#ifndef NDEBUG
	std::cerr << "ok" << std::endl;
#endif

	// *** compute topological hierarchy of polygons
	if (merged_plys) {
		for (std::vector<GEOLIB::Polyline*>::const_iterator it(merged_plys->begin());
			it!=merged_plys->end(); it++) {
			if ((*it)->isClosed()) {
				_polygon_tree_list.push_back(new GMSHPolygonTree(new GEOLIB::Polygon(*(*it), true), NULL, _geo_objs, _gmsh_geo_name, _mesh_density_strategy));
			}
		}
		std::cout << "[GMSHInterface::writeGMSHInputFile] compute topological hierarchy - detected "
						<< _polygon_tree_list.size() << " polygons" << std::endl;
		GEOLIB::createPolygonTrees<FileIO::GMSHPolygonTree>(_polygon_tree_list);
		std::cout << "[GMSHInterface::writeGMSHInputFile] compute topological hierarchy - calculated "
								<< _polygon_tree_list.size() << " polygon trees" << std::endl;
	} else {
		return;
	}

	// *** insert stations and polylines (except polygons) in the appropriate object of
	//     class GMSHPolygonTree
	if (merged_stations) {
		const size_t n_stations(merged_stations->size());
		for (size_t k(0); k < n_stations; k++) {
			bool found(false);
			for (std::list<GMSHPolygonTree*>::iterator it(_polygon_tree_list.begin());
							it != _polygon_tree_list.end() && !found; it++) {
				if ((*it)->insertStation((*merged_stations)[k])) {
					found = true;
				}
			}
		}
	}

	const size_t n_plys(merged_plys->size());
	for (size_t k(0); k<n_plys; k++) {
		if (! (*merged_plys)[k]->isClosed()) {
			for (std::list<GMSHPolygonTree*>::iterator it(_polygon_tree_list.begin());
				it != _polygon_tree_list.end(); it++) {
				(*it)->insertPolyline(new GEOLIB::PolylineWithSegmentMarker(*(*merged_plys)[k]));
			}
		}
	}

	// *** init mesh density strategies
	for (std::list<GMSHPolygonTree*>::iterator it(_polygon_tree_list.begin());
		it != _polygon_tree_list.end(); it++) {
		(*it)->initMeshDensityStrategy();
	}

	// *** create GMSH data structures
	_gmsh_pnts.resize(merged_pnts->size());
	for (std::list<GMSHPolygonTree*>::iterator it(_polygon_tree_list.begin());
		it != _polygon_tree_list.end(); it++) {
		(*it)->createGMSHPoints(_gmsh_pnts);
	}

	// *** finally write data :-)
	writePoints(out);
	for (std::list<GMSHPolygonTree*>::iterator it(_polygon_tree_list.begin());
		it != _polygon_tree_list.end(); it++) {
		(*it)->writeLineLoop(_n_lines, _n_plane_sfc, out);
		(*it)->writeSubPolygonsAsLineConstraints(_n_lines, _n_plane_sfc-1, out);
		(*it)->writeLineConstraints(_n_lines, _n_plane_sfc-1, out);
//		(*it)->writeStations(_n_plane_sfc-1, out);
	}
//	writeLines(out);

	_geo_objs.removeSurfaceVec(_gmsh_geo_name);
	_geo_objs.removePolylineVec(_gmsh_geo_name);
	_geo_objs.removePointVec(_gmsh_geo_name);
	_geo_objs.removeStationVec(_gmsh_station_name);
}

//size_t GMSHInterface::createGMSHPointsInsideBoundingPolygon(
//				GEOLIB::Polygon const* const bounding_polygon,
//				const std::vector<GEOLIB::Point*> &pnt_vec, bool station)
//{
//	// write points
//	const size_t n(pnt_vec.size());
//	size_t cnt(0);
//	if (!station || dynamic_cast<GMSHAdaptiveMeshDensity*> (_mesh_density_strategy) == NULL) {
//		for (size_t k(0); k < n; k++) {
//			if (bounding_polygon->isPntInPolygon(*(pnt_vec[k]))) {
//				_gmsh_pnts.push_back(GMSHPoint((*(pnt_vec[k])), _n_pnt_offset + cnt,
//								_mesh_density_strategy->getMeshDensityAtPoint(pnt_vec[k])));
//				_gmsh_pnt_id_to_geo_pnt_id_map[k] = _n_pnt_offset + cnt;
//				cnt++;
//			} else {
//				_gmsh_pnt_id_to_geo_pnt_id_map[k] = std::numeric_limits<size_t>::max();
//			}
//		}
//	} else {
//		for (size_t k(0); k < n; k++) {
//			if (bounding_polygon->isPntInPolygon(*(pnt_vec[k]))) {
//				_gmsh_pnts.push_back(GMSHPoint((*(pnt_vec[k])), _n_pnt_offset + cnt,
//									dynamic_cast<GMSHAdaptiveMeshDensity*> (_mesh_density_strategy)->getMeshDensityAtStation(pnt_vec[k])));
//				_gmsh_pnt_id_to_geo_pnt_id_map[k] = _n_pnt_offset + cnt;
//				cnt++;
//			} else {
//				_gmsh_pnt_id_to_geo_pnt_id_map[k] = std::numeric_limits<size_t>::max();
//			}
//		}
//	}
//	_n_pnt_offset += cnt;
//	return cnt;
//}

void GMSHInterface::writePoints(std::ostream& out) const
{
	const size_t n_gmsh_pnts(_gmsh_pnts.size());
	for (size_t k(0); k<n_gmsh_pnts; k++) {
		out << *(_gmsh_pnts[k]) << std::endl;
	}
}

//void GMSHInterface::writeLines(std::ostream& out) const
//{
//	const size_t n_gmsh_line_loops(_gmsh_line_loops.size());
//	_gmsh_line_loops[0]->setSurface(true);
//	for (size_t k(0); k<n_gmsh_line_loops; k++) {
//		_gmsh_line_loops[k]->write(out, _n_lines, _n_plane_sfc);
//	}
//}

//void GMSHInterface::writeGMSHPolyline(std::vector<GEOLIB::Polyline*> const& ply_vec, size_t ply_id,
//				size_t offset, std::ostream& out)
//{
//	size_t local_offset(this->_n_pnt_offset - offset);
//	GEOLIB::Polyline const* const ply(ply_vec[ply_id]);
//	size_t s(ply->getNumberOfPoints());
//
//	_ply_id_gmsh_line_mapping.push_back(PolylineGMSHMapping(ply_id, _n_lines, _n_lines + s - 2));
//
//	// write line segments (= Line) of the polyline
//	for (size_t j(0); j < s - 1; j++) {
//		out << "Line(" << _n_lines + j << ") = {" << local_offset
//						+ _gmsh_pnt_id_to_geo_pnt_id_map[ply->getPointID(j)] << "," << local_offset
//						+ _gmsh_pnt_id_to_geo_pnt_id_map[ply->getPointID(j + 1)] << "};"
//						<< std::endl;
//	}
//	// write the line segments contained in the polyline (=Line Loop)
//	out << "Line Loop (" << _n_lines + s - 1 << ") = {";
//	for (size_t j(0); j < s - 2; j++)
//		out << _n_lines + j << ",";
//	out << _n_lines + s - 2 << "};" << std::endl;
//	_n_lines += s;
//}

//void GMSHInterface::writeGMSHPolylines(const std::vector<GEOLIB::Polyline*>& ply_vec,
//				std::ostream &out)
//{
//	size_t n(ply_vec.size());
//	for (size_t k(0); k < n; k++) {
//		// write k-th polyline
//		writeGMSHPolyline(ply_vec, k, 0, out);
//	}
//	out << std::endl;
//}

//size_t GMSHInterface::writeGMSHPolygon(const std::vector<GEOLIB::Polyline*>& ply_vec, size_t id,
//				size_t offset, std::ostream & out)
//{
//	writeGMSHPolyline(ply_vec, id, offset, out);
//	return _n_lines - 1;
//}

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

//std::list<size_t> GMSHInterface::findHolesInsidePolygon(const std::vector<GEOLIB::Polyline*>* plys,
//				size_t i) const
//{
//	GEOLIB::Polygon polygon(*((*plys)[i]));
//	std::list<size_t> polygon_list;
//	polygon_list.push_back(i);
//	for (size_t j = 0; j < plys->size(); j++) {
//		// check if polygons are located completely inside the given polygon
//		if ((i != j) && ((*plys)[j]->isClosed())) {
//			GEOLIB::Polyline* line((*plys)[j]);
//			bool isInside(true);
//			for (size_t k = 0; k < line->getNumberOfPoints(); k++) {
//				if (!polygon.isPntInPolygon(*(line->getPoint(k)))) {
//					isInside = false;
//					break;
//				}
//			}
//			if (isInside) polygon_list.push_back(j);
//		}
//	}
//	return polygon_list;
//}

//bool GMSHInterface::isPolygonInOtherPolygon(const std::vector<GEOLIB::Polyline*>* plys, size_t i, size_t &j) const
//{
//	if (!((*plys)[i])->isClosed()) return false;
//
//	const size_t n_plys(plys->size());
//	bool ply_i_is_in_ply_j(false);
//	for (j = 0; j < n_plys && !ply_i_is_in_ply_j; j++) {
//		if ((i != j) && (((*plys)[j])->isClosed())) {
//			GEOLIB::Polygon const& polygon_j(*((*plys)[j]));
//			ply_i_is_in_ply_j = polygon_j.isPolylineInPolygon(*((*plys)[i]));
//		}
//	}
//
//	if (i == j) return false;
//	return ply_i_is_in_ply_j;
//}

//void GMSHInterface::writePlaneSurface(std::list<size_t> const& polygon_list, bool respect_holes)
//{
//	_out << "Plane Surface (" << _n_plane_sfc << ") = {" << std::flush;
//	std::list<size_t>::const_iterator it(polygon_list.begin());
//	_out << *it << std::flush;
//	if (respect_holes) {
//		for (++it; it != polygon_list.end(); ++it)
//			_out << ", " << *it << std::flush;
//	}
//	_out << "};" << std::endl;
//	_n_plane_sfc++;
//}

//GEOLIB::Polygon* GMSHInterface::getBoundingPolygon(
//				std::vector<GEOLIB::Polyline*> const & all_polylines, size_t& bp_idx) const
//{
//	GEOLIB::Polygon* bounding_polygon(NULL);
//	const size_t n_polylines(all_polylines.size());
//	for (size_t k(0); k < n_polylines; k++) {
//		if (all_polylines[k]->isClosed()) // == Polygon
//		{
//			if (bounding_polygon) { // we have already a bounding polygon
//				if (!bounding_polygon->isPolylineInPolygon(*(all_polylines[k]))) {
//					GEOLIB::Polygon* tmp_polygon(new GEOLIB::Polygon(*(all_polylines[k])));
//					if (tmp_polygon->isPolylineInPolygon(*bounding_polygon)) {
//						// found new bounding polygon
//						delete bounding_polygon;
//						bounding_polygon = tmp_polygon;
//						bp_idx = k;
//					} else {
//						std::cerr << "INFO: there is no inclusion relation between the polygons "
//										<< k << " and " << bp_idx << std::endl;
//					}
//				}
//			} else {
//				bounding_polygon = new GEOLIB::Polygon(*(all_polylines[k]));
//				bp_idx = k;
//			}
//		}
//	}
//
//	if (!bounding_polygon) {
//		std::cerr
//						<< "WARNING: GMSHInterface::writeAllDataToGMSHInputFile: did not found bounding polygon - abort writing"
//						<< std::endl;
//		return NULL;
//	}
//
//	return bounding_polygon;
//}

//void GMSHInterface::createBoundingPolygon(GEOLIB::Polygon const* const bounding_polygon,
//				size_t pnt_offset)
//{
//	std::cout << "\tcreating bounding polygon ... " << std::flush;
//
//	// *** create GMSHLineLoop data structure in order to store the bounding polygon
//	GMSHLineLoop* gmsh_bounding_polygon(new GMSHLineLoop(true));
//	_gmsh_line_loops.push_back(gmsh_bounding_polygon);
//
//	// create line segments (= Line) of the polyline
//	size_t s(bounding_polygon->getNumberOfPoints());
//	for (size_t j(0); j < s - 1; j++) {
//		const size_t start(pnt_offset + _gmsh_pnt_id_to_geo_pnt_id_map[bounding_polygon->getPointID(j)]);
//		const size_t end(pnt_offset + _gmsh_pnt_id_to_geo_pnt_id_map[bounding_polygon->getPointID(j + 1)]);
//		gmsh_bounding_polygon->addLine(new GMSHLine(start, end));
//	}
//
//	std::cout << "ok" << std::endl;
//}

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

//void GMSHInterface::getStationPoints(const GEOLIB::GEOObjects &geo,
//				std::vector<GEOLIB::Point*>& station_points) const
//{
//	std::vector<std::string> stn_names;
//	geo.getStationVectorNames(stn_names);
//
//	for (std::vector<std::string>::const_iterator it(stn_names.begin()); it != stn_names.end(); ++it) {
//		const std::vector<GEOLIB::Point*>* pnts(geo.getPointVec(*it));
//		station_points.insert(station_points.end(), pnts->begin(), pnts->end());
//	}
//}

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

//void GMSHInterface::createPolylinesAsConstraints(GEOLIB::Polygon * bounding_polygon, size_t bp_idx,
//				std::vector<GEOLIB::Polyline*> const& all_polylines, std::ostream & out)
//{
//	const size_t n_polylines(all_polylines.size());
//	for (size_t k(0); k < n_polylines; k++) {
//		if (k != bp_idx) {
//			bool begin_line_pnt_inside_polygon(true);
//			bool end_line_pnt_inside_polygon(true);
//
//			size_t n_pnts_in_ply(all_polylines[k]->getNumberOfPoints());
//
//			// write line segments (= Line) of the polyline
//			for (size_t j(0); j < n_pnts_in_ply - 1; j++) {
//				// check if line segment is contained in bounding polygon
//				bool line_seg_is_already_used(GEOLIB::containsEdge(
//								*(dynamic_cast<GEOLIB::Polyline*> (bounding_polygon)),
//								(all_polylines[k])->getPointID(j), (all_polylines[k])->getPointID(j
//												+ 1)));
//				if (line_seg_is_already_used) {
//					std::cout << "line segment intersects bp" << std::endl;
//				}
//				// check if line segment is contained in a previous polyline
//				for (size_t i(0); i < k && !line_seg_is_already_used; i++)
//					line_seg_is_already_used = GEOLIB::containsEdge(*(all_polylines[i]),
//									(all_polylines[k])->getPointID(j),
//									(all_polylines[k])->getPointID(j + 1));
//
//				if (!line_seg_is_already_used) {
//					// check if first point of line segment is inside bounding polygon
//					if (j == 0) {
//						if (k == 9) std::cout << "pnt in polygon test" << std::endl;
//						begin_line_pnt_inside_polygon = bounding_polygon->isPntInPolygon(
//										*(all_polylines[k])->getPoint(j));
//					}
//					// check if end point of the line is inside bounding polygon
//					end_line_pnt_inside_polygon = bounding_polygon->isPntInPolygon(
//									*(all_polylines[k])->getPoint(j + 1));
//
//					if (begin_line_pnt_inside_polygon && end_line_pnt_inside_polygon) {
//						if (!GEOLIB::isLineSegmentIntersecting(*bounding_polygon,
//										*(all_polylines[k])->getPoint(j),
//										*(all_polylines[k])->getPoint(j + 1))) {
//							out << "Line(" << _n_lines + j << ") = {"
//											<< _gmsh_pnt_id_to_geo_pnt_id_map[(all_polylines[k])->getPointID(
//															j)] << ","
//											<< _gmsh_pnt_id_to_geo_pnt_id_map[(all_polylines[k])->getPointID(
//															j + 1)] << "};" << std::endl;
//							// write line as constraint
//							out << "Line {" << _n_lines + j << "} In Surface {" << _n_plane_sfc - 1
//											<< "};" << std::endl;
//						}
//						//						else {
//						//							std::cout << "segment [(" << *(all_polylines[k])->getPoint(j) <<"), (" <<  *(all_polylines[k])->getPoint(j+1) << ")] intersects with bp" << std::endl;
//						//						}
//					} else {
//						if (begin_line_pnt_inside_polygon && !end_line_pnt_inside_polygon) {
//							// create new point
//							size_t seg_num(0);
//							GEOLIB::Point * s(bounding_polygon->getIntersectionPointPolygonLine(
//											*(all_polylines[k])->getPoint(j),
//											*(all_polylines[k])->getPoint(j + 1), seg_num));
//							if (s != NULL) {
//								std::cout.precision(20);
//								std::cout << "Point(" << _n_pnt_offset << ") = " << *s << std::endl;
//								out << "Point(" << _n_pnt_offset << ") = {" << (*s)[0] << ","
//												<< (*s)[1] << "," << (*s)[2];
//								out << ", " << _mesh_density_strategy->getMeshDensityAtPoint(
//												all_polylines[k]->getPoint(j));
//								out << "}; // new point of polyline " << k << std::endl;
//								// write line
//								out << "Line(" << _n_lines + j << ") = {"
//												<< _gmsh_pnt_id_to_geo_pnt_id_map[(all_polylines[k])-> getPointID(
//																j)] << "," << _n_pnt_offset << "};"
//												<< std::endl;
//								// write line as constraint
//								out << "Line {" << _n_lines + j << "} In Surface {" << _n_plane_sfc
//												- 1 << "};" << std::endl;
//								_n_pnt_offset++;
//								delete s;
//							}
//						}
//						if (!begin_line_pnt_inside_polygon && end_line_pnt_inside_polygon) {
//							// create new point
//							size_t seg_num(0);
//							GEOLIB::Point * s(bounding_polygon->getIntersectionPointPolygonLine(
//											*(all_polylines[k])->getPoint(j),
//											*(all_polylines[k])->getPoint(j + 1), seg_num));
//							if (s != NULL) {
//								out << "Point(" << _n_pnt_offset << ") = {" << (*s)[0] << ","
//												<< (*s)[1] << "," << (*s)[2];
//								out << ", " << _mesh_density_strategy->getMeshDensityAtPoint(
//												all_polylines[k]->getPoint(j));
//								out << "}; // new point of polyline " << k << std::endl;
//								// write line
//								out << "Line(" << _n_lines + j << ") = {" << _n_pnt_offset << ","
//												<< _gmsh_pnt_id_to_geo_pnt_id_map[(all_polylines[k])->getPointID(
//																j + 1)] << "};" << std::endl;
//								// write line as constraint
//								out << "Line {" << _n_lines + j << "} In Surface {" << _n_plane_sfc
//												- 1 << "};" << std::endl;
//								_n_pnt_offset++;
//								delete s;
//							}
//						}
//					}
//					begin_line_pnt_inside_polygon = end_line_pnt_inside_polygon;
//				}
//				//				for (size_t j(0); j < s - 1; j++) {
//				//					// write line as constraint
//				//					out << "Line {" << _n_lines + j << "} In Surface {" << _n_plane_sfc - 1 << "};" << std::endl;
//				//				}
//			}
//			// update line counter
//			_n_lines += n_pnts_in_ply;
//		}
//	}
//}

} // end namespace FileIO
