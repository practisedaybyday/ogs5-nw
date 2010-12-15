/*
 * GMSHInterface.cpp
 *
 *  Created on: Apr 29, 2010
 *      Author: TF
 */

#include <fstream>
#include <vector>
#include <list>

// Base
#include "swap.h"

// FileIO
#include <GMSHInterface.h>

// GEOLIB
#include "Point.h"
#include "Polyline.h"
#include "Polygon.h"
#include "SimplePolygonTree.h"
#include "QuadTree.h"

namespace FileIO {

GMSHInterface::GMSHInterface (const std::string &fname) :
	_n_pnt_offset (0),
	_n_lines (0), // _n_line_loops (0),
	_n_plane_sfc(0)
{
	// open file
	_out.open (fname.c_str());
	// check file stream
	if (!_out) {
		std::cerr << "could not open file " << fname << std::endl;
		return;
	}

	_out << "// GMSH input file created by OpenGeoSys" << std::endl;
	_out << std::endl;
	_out.precision (20);
}

GMSHInterface::~GMSHInterface ()
{
	_out.close ();
}

void GMSHInterface::writeGMSHPoints(const std::vector<GEOLIB::Point*>& pnt_vec)
{
	// write points
//	float characteristic_len (10.0);
	size_t n (pnt_vec.size());
	for (size_t k(0); k<n; k++) {
		_out << "Point(" << _n_pnt_offset + k << ") = {" << (*(pnt_vec[k]))[0] << ","
			<< (*(pnt_vec[k]))[1] << "," << (*(pnt_vec[k]))[2]
//			<< "," << characteristic_len
			<< "};" << std::endl;
	}
	_n_pnt_offset += n;
}

void GMSHInterface::writeGMSHPolyline (const GEOLIB::Polyline* ply)
{
	size_t s (ply->getNumberOfPoints());
	// write line segments (= Line) of the polyline
	for (size_t j(0); j<s-1; j++) {
		_out << "Line(" << _n_lines+j << ") = {" << _n_pnt_offset + ply->getPointID(j) << ","
				<< _n_pnt_offset + ply->getPointID(j+1) << "};" << std::endl;
	}
	// write the line segments contained in the polyline (=Line Loop)
	_out << "Line Loop (" << _n_lines + s - 1 << ") = {";
	for (size_t j(0); j<s-2; j++) {
		_out << _n_lines + j << ",";
	}
	_out << _n_lines + s-2 << "};" << std::endl;
	_n_lines += s;
}

void GMSHInterface::writeGMSHPolylines(const std::vector<GEOLIB::Polyline*>& ply_vec)
{
	size_t n (ply_vec.size());
	for (size_t k(0); k<n; k++) {
		// write k-th polyline
		writeGMSHPolyline (ply_vec[k]);
	}
	_out << std::endl;
}

void GMSHInterface::writeGMSHPolygon(const GEOLIB::Polygon& polygon)
{
	writeGMSHPolyline (&polygon);
	_polygon_list.push_back (_n_lines-1);
}

bool GMSHInterface::writeGMSHInputFile(const std::string &proj_name, const GEOLIB::GEOObjects& geo)
{
	std::cerr << "GMSHInterface::writeGMSHInputFile " << std::endl;
	std::cerr << "get data from geo ... " << std::flush;
	// get data from geo
	const std::vector<GEOLIB::Point*> *pnts (geo.getPointVec (proj_name));
	const std::vector<GEOLIB::Polyline*> *plys (geo.getPolylineVec (proj_name));
	std::cerr << "ok" << std::endl;

	// check file stream
	if (!_out) return false;

	// write points
	writeGMSHPoints (*pnts);

	// write Polylines
	writeGMSHPolylines (*plys);
	std::cerr << "ok" << std::endl;

	return true;
}

void GMSHInterface::writePlaneSurface ()
{
	_out << "Plane Surface (" << _n_plane_sfc << ") = {" << std::flush;
	std::list<size_t>::const_iterator it (_polygon_list.begin());
	_out << *it << std::flush;
	for (it++; it != _polygon_list.end(); it++)
		_out << ", " << *it << std::flush;
	_out << "};" << std::endl;
	_n_plane_sfc++;
	_polygon_list.clear();
}

void GMSHInterface::writeAllDataToGMSHInputFile (GEOLIB::GEOObjects& geo)
{
	// check file stream
	if (! _out) return;

	std::cerr << "GMSHInterface::writeGMSHInputFile " << std::endl;

	// get names of all available data sources
	std::vector<std::string> geo_names;
	geo.getGeometryNames(geo_names);
	// get station names
	std::vector<std::string> geo_station_names;
	geo.getStationNames(geo_station_names);

	size_t pnt_offset (0);
	std::vector<GEOLIB::Point*> all_points;
	std::vector<GEOLIB::Polyline*> all_polylines;
	// fetch points and polylines and add them to the vectors, add points to the QuadTree
	for (std::vector<std::string>::const_iterator it (geo_names.begin());
			it != geo_names.end(); it++) {
		// get data from geo
		std::cout << "fetch data for " << *it << " " << std::flush;
		const std::vector<GEOLIB::Point*> *pnts (geo.getPointVec (*it));
		const std::vector<GEOLIB::Polyline*> *plys (geo.getPolylineVec (*it));
		std::cerr << "ok" << std::endl;

		// insert points into vector all_points
		all_points.insert (all_points.end(), pnts->begin(), pnts->end());

		for (size_t k(0); k<plys->size(); k++) {
			size_t pos (all_polylines.size());
			// insert new polyline
			all_polylines.push_back (new GEOLIB::Polyline (all_points));
			// copy points
			for (size_t j(0); j<(*plys)[k]->getNumberOfPoints(); j++) {
				// set points of polyline
				(all_polylines[pos])->addPoint (pnt_offset + ((*plys)[k])->getPointID(j));
			}
		}
		pnt_offset += pnts->size();
	}

	std::vector<GEOLIB::Point*> all_stations;
	for (std::vector<std::string>::const_iterator it (geo_station_names.begin());
		it != geo_station_names.end(); it++) {
		// get data from geo
		std::cout << "fetch data for " << *it << " " << std::flush;
		const std::vector<GEOLIB::Point*> *pnts (geo.getPointVec (*it));
		std::cerr << "ok" << std::endl;
		// insert points into vector all_points
		all_stations.insert (all_stations.end(), pnts->begin(), pnts->end());
	}

	// *** QuadTree - determining bounding box
	std::cout << "computing axis aligned bounding box for quadtree ... " << std::flush;
	// determine axis aligned bounding box
	GEOLIB::Point ll(std::numeric_limits<double>::max(), std::numeric_limits<double>::max(), 0);
	GEOLIB::Point ur(std::numeric_limits<double>::min(), std::numeric_limits<double>::min(), 0);
	for (size_t k(0); k<all_points.size(); k++) {
		if ((*(all_points[k]))[0] < ll[0]) ll[0] = (*(all_points[k]))[0];
		if ((*(all_points[k]))[1] < ll[1]) ll[1] = (*(all_points[k]))[1];
		if ((*(all_points[k]))[0] > ur[0]) ur[0] = (*(all_points[k]))[0];
		if ((*(all_points[k]))[1] > ur[1]) ur[1] = (*(all_points[k]))[1];
	}
	std::cout << "ok" << std::endl;
	// *** QuadTree - create object
	const size_t number_of_point_per_quadtree_node (2);
	std::cout << "creating quadtree ... " << std::flush;
	GEOLIB::QuadTree<GEOLIB::Point, number_of_point_per_quadtree_node> quad_tree (ll, ur);
	std::cout << "ok" << std::endl;

	// *** QuadTree - insert points
	std::cout << "inserting " << all_points.size() << " points into quadtree ... " << std::flush;
	for (size_t k(0); k < all_points.size(); k++) {
		quad_tree.addPoint (all_points[k]);
	}
	std::cout << "ok" << std::endl;

	// *** QuadTree - insert stations
	std::cout << "inserting " << all_stations.size() << " stations into quadtree ... " << std::flush;
	for (size_t k(0); k<all_stations.size(); k++) {
		quad_tree.addPoint (all_stations[k]);
	}
	std::cout << "ok" << std::endl;

	// *** QuadTree - balance
	std::cout << "balancing quadtree ... " << std::flush;
	quad_tree.balance ();
	std::cout << "ok" << std::endl;

	// *** GMSH - write all non-station points
	const size_t n (all_points.size());
	for (size_t k(0); k<n; k++) {
		GEOLIB::Point ll, ur;
		quad_tree.getLeaf (*(all_points[k]), ll, ur);
		double mesh_density (0.3*(ur[0]-ll[0])); // scaling with 0.3 - do not know if this is a good value
		_out << "Point(" << _n_pnt_offset + k << ") = {" << (*(all_points[k]))[0] << ","
			<< (*(all_points[k]))[1] << "," << (*(all_points[k]))[2]
			<< "," << mesh_density
			<< "};" << std::endl;
	}
	_n_pnt_offset += n;

	// search bounding polygon
	GEOLIB::Polygon* bounding_polygon (NULL);
	size_t bp_idx (0); // bounding polygon index
	const size_t n_polylines (all_polylines.size());
	for (size_t k(0); k<n_polylines; k++) {
		if (all_polylines[k]->isClosed ()) { // == Polygon
			if (bounding_polygon) { // we have already a bounding polygon
				if (! bounding_polygon->isPolylineInPolygon (*(all_polylines[k]))) {
					GEOLIB::Polygon* tmp_polygon (new GEOLIB::Polygon (*(all_polylines[k])));
					if (tmp_polygon->isPolylineInPolygon (*bounding_polygon)) {
						// found new bounding polygon
						delete bounding_polygon;
						bounding_polygon = tmp_polygon;
						bp_idx = k;
					} else {
						std::cout << "there is no inclusion relation between the polygons " << k << " and " << bp_idx << std::endl;
					}
				}
			} else {
				bounding_polygon = new GEOLIB::Polygon (*(all_polylines[k]));
				bp_idx = k;
			}
		}
	}
//	std::cout << "bounding polygon is polygon " << bp_idx << std::endl;

	std::cout << "write bounding polygon ... " << std::flush;
	// write bounding polygon
	size_t s (bounding_polygon->getNumberOfPoints());
	// write line segments (= Line) of the polyline
	for (size_t j(0); j<s-1; j++) {
		_out << "Line(" << _n_lines+j << ") = {" <<  bounding_polygon->getPointID(j) << ","
				<< bounding_polygon->getPointID(j+1) << "};" << std::endl;
	}
	// write the line segments contained in the polyline (=Line Loop)
	_out << "Line Loop (" << _n_lines + s - 1 << ") = {";
	for (size_t j(0); j<s-2; j++) {
		_out << _n_lines + j << ",";
	}
	_out << _n_lines + s-2 << "};" << std::endl;
	_n_lines += s;
	// write plane surface
	_out << "Plane Surface (" << _n_plane_sfc << ") = {" << _n_lines-1 << "};" << std::endl;
	_n_plane_sfc++;
	std::cout << "ok" << std::endl;

	// write all other polylines as constraints
	for (size_t k(0); k<n_polylines; k++) {
		if (k != bp_idx) {
			size_t s (all_polylines[k]->getNumberOfPoints());

			// write line segments (= Line) of the polyline
			for (size_t j(0); j<s-1; j++) {
				// check if line segment is contained in bounding polygon
				bool line_seg_is_already_used (GEOLIB::containsEdge (*(dynamic_cast<GEOLIB::Polyline*>(bounding_polygon)), (all_polylines[k])->getPointID(j), (all_polylines[k])->getPointID(j+1)));
				// check if line segment is contained in a previous polyline
				for (size_t i(0); i<k && !line_seg_is_already_used; i++) {
					line_seg_is_already_used = GEOLIB::containsEdge (*(all_polylines[i]), (all_polylines[k])->getPointID(j), (all_polylines[k])->getPointID(j+1));
				}

				if (!line_seg_is_already_used) {
					_out << "Line(" << _n_lines+j << ") = {" << (all_polylines[k])->getPointID(j) << ","
							<< (all_polylines[k])->getPointID(j+1) << "};" << std::endl;
					// write line as constraint
					_out << "Line {" << _n_lines + j << "} In Surface {" << _n_plane_sfc-1 << "};" << std::endl;
				}
			}
			// update line counter
			_n_lines += s;
		}
	}

	// write stations as constraints
	_out << "// Stations" << std::endl;
	const size_t n_stations (all_stations.size());
	for (size_t k(0); k<n_stations; k++) {
		GEOLIB::Point ll, ur;
		quad_tree.getLeaf (*(all_stations[k]), ll, ur);
		double mesh_density (0.3*(ur[0]-ll[0])); // scaling with 0.3 - do not know if this is a good value
		_out << "Point(" << _n_pnt_offset + k << ") = {" << (*(all_stations[k]))[0] << ","
			<< (*(all_stations[k]))[1] << "," << (*(all_stations[k]))[2] << "," << mesh_density
			<< "};" << std::endl;

		_out << "Point {" << pnt_offset+k << "} In Surface {" << _n_plane_sfc-1 << "};" << std::endl;
	}
	_n_pnt_offset += n_stations;

	// write Steiner points
	std::list<GEOLIB::QuadTree<GEOLIB::Point, number_of_point_per_quadtree_node>*> leaf_list;
	quad_tree.getLeafs (leaf_list);
	_out << "// write Steiner points" << std::endl;
	for (std::list<GEOLIB::QuadTree<GEOLIB::Point, number_of_point_per_quadtree_node>*>::const_iterator it (leaf_list.begin());
		it != leaf_list.end(); it++) {
		if ((*it)->getPoints().empty()) {
			// compute point from square
			GEOLIB::Point ll, rr;
			(*it)->getSquarePoints (ll, rr);
			GEOLIB::Point mid_point (0.5*(rr[0]+ll[0]), 0.5*(rr[1]+ll[1]), 0.5*(rr[2]+ll[2]));
//			std::cout << "insert Steiner point: " << mid_point << std::endl;
			if (bounding_polygon->isPntInPolygon (mid_point)) {
				_out << "Point(" << _n_pnt_offset << ") = {" << mid_point[0] << ","
							<< mid_point[1] << "," << mid_point[2]
							<< "," << 0.5*(rr[0]-ll[0])
							<< "};" << std::endl;
				_out << "Point {" << _n_pnt_offset << "} In Surface {" << _n_plane_sfc-1 << "};" << std::endl;
				_n_pnt_offset++;
			}

		}
	}
	std::cout << "ok" << std::endl;

/*
	// debug
	std::vector<GEOLIB::Point*> *quad_pnts (new std::vector<GEOLIB::Point*>);
	std::vector<GEOLIB::Polyline*> *quad_plys (new std::vector<GEOLIB::Polyline*>);
	quad_tree.getQuadTree (*quad_pnts, *quad_plys);
	std::string name ("QuadTree");
	geo.addPointVec (quad_pnts, name, NULL);
	geo.addPolylineVec (quad_plys, name);
*/

	//	// go through all geometric data sets and get polygons
	//	std::vector<GEOLIB::Polygon*> polygon_vec;
	//	for (std::vector<std::string>::const_iterator it (geo_names.begin());
	//		it != geo_names.end(); it++) {
	//		const std::vector<GEOLIB::Polyline*> *polylines (geo.getPolylineVec (*it));
	//		for (std::vector<GEOLIB::Polyline*>::const_iterator it (polylines->begin()); it != polylines->end(); it++) {
	//			if ((*it)->isClosed ()) {
	//				polygon_vec.push_back (new GEOLIB::Polygon (*(*it)));
	//			}
	//		}
	//	}
	//
	//	// we assume that all polygons are simple polygons
	//	// forest consist of (hierarchy) trees
	//	std::list<GEOLIB::SimplePolygonTree*> polygon_forest;
	//	// create polygon forest
	//	for (std::vector<GEOLIB::Polygon*>::iterator
	//			polygon_it(polygon_vec.begin()); polygon_it != polygon_vec.end(); polygon_it++) {
	//		// get the list and insert the elements as SimplePolygonTree items into the forest
	//		const std::list<GEOLIB::Polygon*> simple_polygon_list(
	//				(*polygon_it)->getListOfSimplePolygons());
	//		for (std::list<GEOLIB::Polygon*>::const_iterator simple_polygon_it(
	//				simple_polygon_list.begin()); simple_polygon_it
	//				!= simple_polygon_list.end(); simple_polygon_it++) {
	//			GEOLIB::SimplePolygonTree *spt(new GEOLIB::SimplePolygonTree(*simple_polygon_it));
	//			polygon_forest.push_back(spt);
	//		}
	//	}
	//
	//	// create the hierarchy
	//	GEOLIB::createPolygonTree(polygon_forest);
	//	std::cout << "\"Polygon forest\" consists of "
	//			<< polygon_forest.size() << " trees" << std::endl;
	//
	//	// *** insert additional Points (for instance Stations) and Polylines (Wadis, rivers, ...)
	//	// Stations
	//	for (std::vector<std::string>::const_iterator it (geo_names.begin());
	//		it != geo_names.end(); it++) {
	//		const std::vector<GEOLIB::Point*> *stations (geo.getStationVec (*it));
	//		if (stations) {
	//			// go through all Stations
	//			for (std::vector<GEOLIB::Point*>::const_iterator it (stations->begin()); it != stations->end(); it++) {
	//				bool nfound (true);
	//				// go through all top level polygons / SimplePolygonTrees
	//				for (std::list<GEOLIB::SimplePolygonTree*>::iterator polygon_it (polygon_forest.begin());
	//					polygon_it != polygon_forest.end() && nfound; polygon_it++) {
	//					if ((*polygon_it)->isGeoObjInside (*it)) {
	//						(*polygon_it)->insertGeoObj (*it);
	//						nfound = false;
	//					}
	//				}
	//			}
	//		}
	//	}
	//	// Polylines
	//	for (std::vector<std::string>::const_iterator it (geo_names.begin());
	//		it != geo_names.end(); it++) {
	//		const std::vector<GEOLIB::Polyline*> *polylines (geo.getPolylineVec (*it));
	//		for (std::vector<GEOLIB::Polyline*>::const_iterator it (polylines->begin()); it != polylines->end(); it++) {
	//			if (! (*it)->isClosed ()) {
	//				bool nfound (true);
	//				// go through all top level polygons / SimplePolygonTrees
	//				for (std::list<GEOLIB::SimplePolygonTree*>::iterator polygon_it (polygon_forest.begin());
	//					polygon_it != polygon_forest.end() && nfound; polygon_it++) {
	//					if ((*polygon_it)->isGeoObjInside (*it)) {
	//						(*polygon_it)->insertGeoObj (*it);
	//						nfound = false;
	//					}
	//				}
	//			}
	//		}
	//	}

}

bool GMSHInterface::isGMSHMeshFile (const std::string& fname)
{
	std::ifstream input (fname.c_str());

	if (!input) {
		std::cerr << "GMSHInterface::isGMSHMeshFile could not open file " << fname << std::endl;
		return false;
	}

	std::string header_first_line;
	input >> header_first_line;
	if (header_first_line.find ("$MeshFormat") != std::string::npos) {
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

} // end namespace FileIO
