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

// FileIO
#include "GMSHInterface.h"

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
GMSHInterface::GMSHInterface (const std::string &fname) :
	_n_pnt_offset (0),
	_n_lines (0),
	_n_plane_sfc(0)
{
	// open file
	_out.open (fname.c_str());
	// check file stream
	if (!_out)
	{
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

void GMSHInterface::writeGMSHPoints(const std::vector<GEOLIB::Point*> &pnt_vec,
                                    GEOLIB::QuadTree<GEOLIB::Point>* quad_tree)
{
	// write points
//	float characteristic_len (10.0);
	const size_t n (pnt_vec.size());
	for (size_t k(0); k < n; k++)
	{
		_out << "Point(" << _n_pnt_offset + k << ") = {" << (*(pnt_vec[k]))[0] << ","
		     << (*(pnt_vec[k]))[1] << "," << (*(pnt_vec[k]))[2];

		if (quad_tree)
		{
			GEOLIB::Point ll, ur;
			quad_tree->getLeaf(*(pnt_vec[k]), ll, ur);
			_out << "," << (0.5 * (ur[0] - ll[0]));
		}

		_out << "};" << std::endl;
	}
	_n_pnt_offset += n;
}

void GMSHInterface::writeGMSHPolyline (const GEOLIB::Polyline* ply, const size_t offset)
{
	size_t local_offset (this->_n_pnt_offset - offset);
	size_t s (ply->getNumberOfPoints());
	// write line segments (= Line) of the polyline
	for (size_t j(0); j < s - 1; j++)
		_out << "Line(" << _n_lines + j << ") = {" << local_offset + ply->getPointID(j) <<
		","
		     << local_offset + ply->getPointID(j + 1) << "};" << std::endl;
	// write the line segments contained in the polyline (=Line Loop)
	_out << "Line Loop (" << _n_lines + s - 1 << ") = {";
	for (size_t j(0); j < s - 2; j++)
		_out << _n_lines + j << ",";
	_out << _n_lines + s - 2 << "};" << std::endl;
	_n_lines += s;
}

void GMSHInterface::writeGMSHPolylines(const std::vector<GEOLIB::Polyline*>& ply_vec)
{
	size_t n (ply_vec.size());
	for (size_t k(0); k < n; k++)
		// write k-th polyline
		writeGMSHPolyline (ply_vec[k], 0);
	_out << std::endl;
}

size_t GMSHInterface::writeGMSHPolygon(const GEOLIB::Polygon& polygon, const size_t offset)
{
	writeGMSHPolyline (&polygon, offset);
	//_polygon_list.push_back (_n_lines-1);
	return _n_lines - 1;
}

bool GMSHInterface::writeGMSHInputFile(const std::string &proj_name,
                                       const GEOLIB::GEOObjects& geo,
                                       bool useStationsAsContraints)
{
	std::cerr << "GMSHInterface::writeGMSHInputFile " << std::endl;
	std::cerr << "get data from geo ... " << std::flush;
	// get data from geo
	const std::vector<GEOLIB::Point*>* pnts (geo.getPointVec (proj_name));
	const std::vector<GEOLIB::Polyline*>* plys (geo.getPolylineVec (proj_name));
	std::cerr << "ok" << std::endl;

	// check file stream
	if (!_out)
		return false;

	// create a quad tree for generate steiner points. this makes sure that the resulting geometry will be suitable for a FEM mesh.
	GEOLIB::QuadTree<GEOLIB::Point>* quad_tree = this->createQuadTreeFromPoints(*pnts, 2);

	std::vector<GEOLIB::Point*> station_points;
	if (useStationsAsContraints)
	{
		station_points = this->getStationPoints(geo);
		for (size_t k(0); k < station_points.size(); k++)
			quad_tree->addPoint (station_points[k]);
	}

	quad_tree->balance();

	// write points
	writeGMSHPoints (*pnts);

	// write Polylines
	std::map<size_t,size_t> geo2gmsh_polygon_id_map;
	std::map<size_t,size_t> geo2gmsh_surface_id_map;

	for (size_t i = 0; i < plys->size(); i++)
	{
		if ((*plys)[i]->isClosed())
		{
			GEOLIB::Polygon polygon(*((*plys)[i]));
			size_t polygon_id = this->writeGMSHPolygon(polygon, pnts->size());
			geo2gmsh_polygon_id_map[i] = polygon_id;
		}
		else
			this->writeGMSHPolyline((*plys)[i], pnts->size());
	}

	for (size_t i = 0; i < plys->size(); i++)
		if ((*plys)[i]->isClosed())
		{
			std::list<size_t> polygon_list = findHolesInsidePolygon(
			        plys,
			        i,
			        geo2gmsh_polygon_id_map);
			this->writePlaneSurface(polygon_list);
			geo2gmsh_surface_id_map[i] = _n_plane_sfc - 1;
		}

	if (useStationsAsContraints)
		this->addPointsAsConstraints(station_points, *plys, geo2gmsh_surface_id_map);

	std::vector<GEOLIB::Point*> steiner_points = this->getSteinerPoints(quad_tree);
	this->addPointsAsConstraints(steiner_points, *plys, geo2gmsh_surface_id_map);

	delete quad_tree;
	for (size_t i = 0; i < steiner_points.size(); i++)
		delete steiner_points[i];

	std::cerr << "ok" << std::endl;

	return true;
}

std::list<size_t> GMSHInterface::findHolesInsidePolygon(const std::vector<GEOLIB::Polyline*>* plys,
                                                        size_t i,
                                                        std::map<size_t,
                                                                 size_t> geo2gmsh_polygon_id_map)
{
	GEOLIB::Polygon polygon(*((*plys)[i]));
	std::list<size_t> polygon_list;
	polygon_list.push_back(geo2gmsh_polygon_id_map[i]);
	for (size_t j = 0; j < plys->size(); j++) // check if polygons are located completely inside the given polygon

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
				polygon_list.push_back(geo2gmsh_polygon_id_map[j]);
		}
	return polygon_list;
}

void GMSHInterface::writePlaneSurface (std::list<size_t> const& polygon_list)
{
	_out << "Plane Surface (" << _n_plane_sfc << ") = {" << std::flush;
	std::list<size_t>::const_iterator it (polygon_list.begin());
	_out << *it << std::flush;
	for (++it; it != polygon_list.end(); ++it)
		_out << ", " << *it << std::flush;
	_out << "};" << std::endl;
	_n_plane_sfc++;
	//_polygon_list.clear();
}

GEOLIB::QuadTree<GEOLIB::Point>* GMSHInterface::createQuadTreeFromPoints(
        std::vector<GEOLIB::Point*> all_points,
        size_t number_of_point_per_quadtree_node)
{
	// *** QuadTree - determining bounding box
#ifndef NDEBUG
	std::cout << "computing axis aligned bounding box for quadtree ... " << std::flush;
#endif
	// determine axis aligned bounding box
	GEOLIB::Point ll(std::numeric_limits<double>::max(), std::numeric_limits<double>::max(), 0);
	GEOLIB::Point ur(std::numeric_limits<double>::min(), std::numeric_limits<double>::min(), 0);
	for (size_t k(0); k < all_points.size(); k++)
	{
		if ((*(all_points[k]))[0] < ll[0])
			ll[0] = (*(all_points[k]))[0];
		if ((*(all_points[k]))[1] < ll[1])
			ll[1] = (*(all_points[k]))[1];
		if ((*(all_points[k]))[0] > ur[0])
			ur[0] = (*(all_points[k]))[0];
		if ((*(all_points[k]))[1] > ur[1])
			ur[1] = (*(all_points[k]))[1];
	}
#ifndef NDEBUG
	std::cout << "ok" << std::endl;
#endif
	// *** QuadTree - create object
#ifndef NDEBUG
	std::cout << "creating quadtree ... " << std::flush;
#endif
	GEOLIB::QuadTree<GEOLIB::Point>* quad_tree = new GEOLIB::QuadTree<GEOLIB::Point>(
	        ll,
	        ur,
	        number_of_point_per_quadtree_node);
	std::cout << "ok" << std::endl;

	// *** QuadTree - insert points
#ifndef NDEBUG
	std::cout << "inserting " << all_points.size() << " points into quadtree ... " <<
	std::flush;
#endif
	for (size_t k(0); k < all_points.size(); k++)
		quad_tree->addPoint (all_points[k]);

#ifndef NDEBUG
	std::cout << "ok" << std::endl;
#endif

	return quad_tree;
}

void GMSHInterface::writeAllDataToGMSHInputFile (
        GEOLIB::GEOObjects& geo,
        std::vector<std::string> const &
        selected_geometries,
        size_t number_of_point_per_quadtree_node,
        double mesh_density_scaling,
        double mesh_density_scaling_station_pnts)
{
	// check file stream
	if (!_out)
		return;

	std::cout << "GMSHInterface::writeGMSHInputFile adaptive " << std::endl;

	std::vector<GEOLIB::Point*> all_points;
	std::vector<GEOLIB::Polyline*> all_polylines;
	std::vector<GEOLIB::Point*> all_stations;
	fetchGeometries (geo, selected_geometries, all_points, all_polylines, all_stations);

	// search bounding polygon
	size_t bp_idx (0); // bounding polygon index
	GEOLIB::Polygon* bounding_polygon (getBoundingPolygon(all_polylines, bp_idx));
	if (!bounding_polygon)
		return;

	GEOLIB::QuadTree<GEOLIB::Point>* quad_tree = this->createQuadTreeFromPoints(
	        all_points,
	        number_of_point_per_quadtree_node);

	// *** QuadTree - insert stations
#ifndef NDEBUG
	std::cout << "inserting " << all_stations.size() << " stations into quadtree ... " <<
	std::flush;
#endif
	for (size_t k(0); k < all_stations.size(); k++)
		quad_tree->addPoint (all_stations[k]);

#ifndef NDEBUG
	std::cout << "ok" << std::endl;
#endif

	// *** QuadTree - balance
#ifndef NDEBUG
	std::cout << "balancing quadtree ... " << std::flush;
#endif
	quad_tree->balance ();
#ifndef NDEBUG
	std::cout << "ok" << std::endl;
#endif

	// *** GMSH - write all non-station points
	const size_t n (all_points.size());
	for (size_t k(0); k < n; k++)
		if (bounding_polygon->isPntInPolygon (*(all_points[k])))
		{
			GEOLIB::Point ll, ur;
			quad_tree->getLeaf (*(all_points[k]), ll, ur);
			double mesh_density (mesh_density_scaling * (ur[0] - ll[0]));
			_out << "Point(" << _n_pnt_offset + k << ") = {" <<
			(*(all_points[k]))[0] << ","
			     << (*(all_points[k]))[1] << "," << (*(all_points[k]))[2]
			     << "," << mesh_density
			     << "};" << std::endl;
		}

	// write the bounding polygon
	writeBoundingPolygon ( bounding_polygon );

	_n_pnt_offset += n;

	// write all other polylines as constraints
	const size_t n_polylines (all_polylines.size());
	for (size_t k(0); k < n_polylines; k++)
		if (k != bp_idx)
		{
			bool begin_line_pnt_inside_polygon (true);
			bool end_line_pnt_inside_polygon (true);

			size_t s (all_polylines[k]->getNumberOfPoints());

			// write line segments (= Line) of the polyline
			for (size_t j(0); j < s - 1; j++)
			{
//				_out << "Line(" << _n_lines+j << ") = {" << (all_polylines[k])->getPointID(j) << ","
//						<< (all_polylines[k])->getPointID(j+1) << "};" << std::endl;
//				// write line as constraint
//				_out << "Line {" << _n_lines + j << "} In Surface {" << _n_plane_sfc-1 << "};" << std::endl;

				// check if line segment is contained in bounding polygon
				bool line_seg_is_already_used (GEOLIB::containsEdge (
				                                       *(dynamic_cast<GEOLIB::
				                                                      Polyline*>(
				                                                 bounding_polygon)),
				                                       (all_polylines[k])->
				                                       getPointID(
				                                               j),
				                                       (all_polylines[k])->
				                                       getPointID(j + 1)));
				// check if line segment is contained in a previous polyline
				for (size_t i(0); i < k && !line_seg_is_already_used; i++)
					line_seg_is_already_used =
					        GEOLIB::containsEdge (
					                *(all_polylines[i]),
					                (all_polylines[k])->
					                getPointID(
					                        j),
					                (all_polylines[k])->getPointID(j + 1));

				if (!line_seg_is_already_used)
				{
					// check if first point of polyline is inside bounding polygon
					if (j == 0)
						begin_line_pnt_inside_polygon =
						        bounding_polygon->isPntInPolygon (*(
						                                                  all_polylines
						                                                  [
						                                                          k])->getPoint(j));
					// check if end point of the line is inside bounding polygon
					end_line_pnt_inside_polygon =
					        bounding_polygon->isPntInPolygon (
					                *(all_polylines[k])->getPoint(j + 1));

					if (begin_line_pnt_inside_polygon &&
					    end_line_pnt_inside_polygon)
					{
						_out << "Line(" << _n_lines + j << ") = {" <<
						(all_polylines[k])->getPointID(j) << ","
						     << (all_polylines[k])->getPointID(j +
						                                  1) << "};" <<
						std::endl;
						// write line as constraint
						_out << "Line {" << _n_lines + j <<
						"} In Surface {" << _n_plane_sfc - 1 <<
						"};" << std::endl;
					}
					else
					{
						if (begin_line_pnt_inside_polygon &&
						    !end_line_pnt_inside_polygon)
						{
							// create new point
							GEOLIB::Point* s (
							        bounding_polygon->
							        getIntersectionPointPolygonLine (
							                *(all_polylines[k])->
							                getPoint(
							                        j),
							                *(all_polylines[k])->
							                getPoint(j + 1)));
							if (s != NULL)
							{
								// write new point as gmsh geo point with mesh density from existing point
								GEOLIB::Point ll, ur;
								quad_tree->getLeaf (
								        *(all_polylines
								          [k])->
								        getPoint(j),
								        ll,
								        ur);
								double mesh_density (
								        0.3 * (ur[0] - ll[0])); // scaling with 0.3 - do not know if this is a good value
								_out << "Point(" <<
								_n_pnt_offset << ") = {" <<
								(*s)[0] << ","
								<< (*s)[1] << "," <<
								(*s)[2] << "," <<
								mesh_density
								<<
								"}; // new end point of polyline "
								<< k << std::endl;
								// write line
								_out << "Line(" << _n_lines + j <<
								") = {" <<
								(all_polylines[k])->
								getPointID(j) << ","
								              << _n_pnt_offset <<
								"};" <<
								std::endl;
								// write line as constraint
								_out << "Line {" << _n_lines + j <<
								"} In Surface {" <<
								_n_plane_sfc - 1 << "};" <<
								std::endl;
								_n_pnt_offset++;
								delete s;
							}
						}
						if (!begin_line_pnt_inside_polygon &&
						    end_line_pnt_inside_polygon)
						{
							// create new point
							GEOLIB::Point* s (
							        bounding_polygon->
							        getIntersectionPointPolygonLine(
							                *(all_polylines[k])->
							                getPoint(
							                        j),
							                *(all_polylines[k])->
							                getPoint(j + 1)));
							if (s != NULL)
							{
								// write new point as gmsh geo point with mesh density from existing point
								GEOLIB::Point ll, ur;
								quad_tree->getLeaf (
								        *(all_polylines
								          [k])->
								        getPoint(
								                j + 1),
								        ll,
								        ur);
								double mesh_density (
								        0.3 * (ur[0] - ll[0])); // scaling with 0.3 - do not know if this is a good value
								_out << "Point(" <<
								_n_pnt_offset << ") = {" <<
								(*s)[0] << ","
								<< (*s)[1] << "," <<
								(*s)[2] << "," <<
								mesh_density
								<<
								"}; // new end point of polyline "
								<< k << std::endl;
								// write line
								_out << "Line(" << _n_lines + j <<
								") = {" << _n_pnt_offset <<
								"," <<
								(all_polylines[k])->getPointID(
								        j + 1)
								<< "};" << std::endl;
								// write line as constraint
								_out << "Line {" << _n_lines + j <<
								"} In Surface {" <<
								_n_plane_sfc - 1 << "};" <<
								std::endl;
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

	// write stations as constraints
	_out << "// Stations" << std::endl;
	const size_t n_stations (all_stations.size());
	for (size_t k(0); k < n_stations; k++)
		if (bounding_polygon->isPntInPolygon (*(all_stations[k])))
		{
			GEOLIB::Point ll, ur;
			quad_tree->getLeaf (*(all_stations[k]), ll, ur);
			double mesh_density (mesh_density_scaling_station_pnts * (ur[0] - ll[0]));
			_out << "Point(" << _n_pnt_offset + k << ") = {" <<
			(*(all_stations[k]))[0] << ","
			     << (*(all_stations[k]))[1] << "," << (*(all_stations[k]))[2] << "," <<
			mesh_density
			     << "};" << std::endl;
			_out << "Point {" << _n_pnt_offset + k << "} In Surface {" <<
			_n_plane_sfc - 1 << "};" << std::endl;
		}
	_n_pnt_offset += n_stations;

	// write Steiner points
	std::list<GEOLIB::QuadTree<GEOLIB::Point>*> leaf_list;
	quad_tree->getLeafs (leaf_list);
	_out << "// Steiner points" << std::endl;
	for (std::list<GEOLIB::QuadTree<GEOLIB::Point>*>::const_iterator it (leaf_list.begin());
	     it != leaf_list.end(); it++)
		if ((*it)->getPoints().empty())
		{
			// compute point from square
			GEOLIB::Point ll, rr;
			(*it)->getSquarePoints (ll, rr);
			GEOLIB::Point mid_point (0.5 * (rr[0] + ll[0]),
			                         0.5 * (rr[1] + ll[1]),
			                         0.5 * (rr[2] + ll[2]));
			if (bounding_polygon->isPntInPolygon (mid_point))
			{
				std::vector<GEOLIB::Point*> const& leaf_pnts ((*it)->getPoints());
				double z_average (0);
				const size_t n_leaf_pnts (leaf_pnts.size());
				if (n_leaf_pnts > 0)
				{
					for (size_t k(0); k < n_leaf_pnts; k++)
						z_average += (*(leaf_pnts[k]))[2];
					z_average /= n_leaf_pnts;
				}
				else
				{
					GEOLIB::QuadTree<GEOLIB::Point> const* const father (
					        (*it)->getFather());
					if (father)
					{
						std::vector<GEOLIB::Point*> const& pnts_father (
						        father->getPoints());

						if (pnts_father.size() > 0)
						{
							for (size_t k(0); k < pnts_father.size();
							     k++)
								z_average += (*(pnts_father[k]))[2];
							z_average /= pnts_father.size();
							std::cout <<
							"DEBUG: pnts_father > 0 in GMSHInterface::writeAllDataToGMSHInputFile while writing Steiner points"
							          << std::endl;
						}
						else
						{
							GEOLIB::QuadTree<GEOLIB::Point> const*
							const child_ne (father->getChild (
							                        GEOLIB::QuadTree<
							                                GEOLIB::
							                                Point>::NE));
							std::vector<GEOLIB::Point*> const&
							leaf_pnts_ne (child_ne->getPoints());
							for (size_t k(0); k < leaf_pnts_ne.size();
							     k++)
								z_average +=
								        (*(leaf_pnts_ne[k]))[2];

							GEOLIB::QuadTree<GEOLIB::Point> const*
							const child_nw (father->getChild (
							                        GEOLIB::QuadTree<
							                                GEOLIB::
							                                Point>::NW));
							std::vector<GEOLIB::Point*> const&
							leaf_pnts_nw (child_nw->getPoints());
							for (size_t k(0); k < leaf_pnts_nw.size();
							     k++)
								z_average +=
								        (*(leaf_pnts_nw[k]))[2];

							GEOLIB::QuadTree<GEOLIB::Point> const*
							const child_sw (father->getChild (
							                        GEOLIB::QuadTree<
							                                GEOLIB::
							                                Point>::SW));
							std::vector<GEOLIB::Point*> const&
							leaf_pnts_sw (child_sw->getPoints());
							for (size_t k(0); k < leaf_pnts_sw.size();
							     k++)
								z_average +=
								        (*(leaf_pnts_sw[k]))[2];

							GEOLIB::QuadTree<GEOLIB::Point> const*
							const child_se (father->getChild (
							                        GEOLIB::QuadTree<
							                                GEOLIB::
							                                Point>::SE));
							std::vector<GEOLIB::Point*> const&
							leaf_pnts_se (child_se->getPoints());
							for (size_t k(0); k < leaf_pnts_se.size();
							     k++)
								z_average +=
								        (*(leaf_pnts_se[k]))[2];
							size_t n_pnts (
							        leaf_pnts_ne.size() +
							        leaf_pnts_nw.size() +
							        leaf_pnts_sw.size() +
							        leaf_pnts_se.size());
							if (n_pnts > 0)
								z_average /= n_pnts;
							else
								std::cout <<
								"DEBUG: n_pnts == 0 in GMSHInterface::writeAllDataToGMSHInputFile while writing Steiner points"
								          << std::endl;
						}
					}
				}
				if (z_average < -1)
				{
					mid_point[2] = z_average;
					_out << "Point(" << _n_pnt_offset << ") = {" <<
					mid_point[0] << ","
					     << mid_point[1] << "," << mid_point[2]
					     << "," << 0.5 * (rr[0] - ll[0])
					     << "};" << std::endl;
					_out << "Point {" << _n_pnt_offset << "} In Surface {" <<
					_n_plane_sfc - 1 << "};" << std::endl;
					_n_pnt_offset++;
				}
			}
		}
	delete quad_tree;
	std::cout << "ok" << std::endl;
}

void GMSHInterface::writeAllDataToGMSHInputFile (
        GEOLIB::GEOObjects& geo,
        std::vector<std::string> const &
        selected_geometries,
        double mesh_density)
{
	// check file stream
	if (!_out)
		return;

	std::cout << "GMSHInterface::writeGMSHInputFile non adaptive" << std::endl;

	std::vector<GEOLIB::Point*> all_points;
	std::vector<GEOLIB::Polyline*> all_polylines;
	std::vector<GEOLIB::Point*> all_stations;
	fetchGeometries (geo, selected_geometries, all_points, all_polylines, all_stations);

	// search bounding polygon
	size_t bp_idx (0); // bounding polygon index
	GEOLIB::Polygon* bounding_polygon (getBoundingPolygon(all_polylines, bp_idx));
	if (!bounding_polygon)
		return;

	// *** GMSH - write all non-station points
	const size_t n (all_points.size());
	for (size_t k(0); k < n; k++)
		if (bounding_polygon->isPntInPolygon (*(all_points[k])))
			_out << "Point(" << _n_pnt_offset + k << ") = {" <<
			(*(all_points[k]))[0] << ","
			     << (*(all_points[k]))[1] << "," << (*(all_points[k]))[2]
			     << "," << mesh_density
			     << "};" << std::endl;

	std::cout << "write bounding polygon ... " << std::flush;
	// write bounding polygon
	writeBoundingPolygon (bounding_polygon);
	_n_pnt_offset += n;

	// write all other polylines as constraints
	const size_t n_polylines (all_polylines.size());
	for (size_t k(0); k < n_polylines; k++)
		if (k != bp_idx)
		{
			bool begin_line_pnt_inside_polygon (true);
			bool end_line_pnt_inside_polygon (true);

			size_t s (all_polylines[k]->getNumberOfPoints());

			// write line segments (= Line) of the polyline
			for (size_t j(0); j < s - 1; j++)
			{
				// check if line segment is contained in bounding polygon
				bool line_seg_is_already_used (GEOLIB::containsEdge (
				                                       *(dynamic_cast<GEOLIB::
				                                                      Polyline*>(
				                                                 bounding_polygon)),
				                                       (all_polylines[k])->
				                                       getPointID(
				                                               j),
				                                       (all_polylines[k])->
				                                       getPointID(j + 1)));
				// check if line segment is contained in a previous polyline
				for (size_t i(0); i < k && !line_seg_is_already_used; i++)
					line_seg_is_already_used =
					        GEOLIB::containsEdge (
					                *(all_polylines[i]),
					                (all_polylines[k])->
					                getPointID(
					                        j),
					                (all_polylines[k])->getPointID(j + 1));

				if (!line_seg_is_already_used)
				{
					// check if first point of polyline is inside bounding polygon
					if (j == 0)
						begin_line_pnt_inside_polygon =
						        bounding_polygon->isPntInPolygon (*(
						                                                  all_polylines
						                                                  [
						                                                          k])->getPoint(j));
					// check if end point of the line is inside bounding polygon
					end_line_pnt_inside_polygon =
					        bounding_polygon->isPntInPolygon (
					                *(all_polylines[k])->getPoint(j + 1));

					if (begin_line_pnt_inside_polygon &&
					    end_line_pnt_inside_polygon)
					{
						_out << "Line(" << _n_lines + j << ") = {" <<
						(all_polylines[k])->getPointID(j) << ","
						     << (all_polylines[k])->getPointID(j +
						                                  1) << "};" <<
						std::endl;
						// write line as constraint
						_out << "Line {" << _n_lines + j <<
						"} In Surface {" << _n_plane_sfc - 1 <<
						"};" << std::endl;
					}
					else
					{
						if (begin_line_pnt_inside_polygon &&
						    !end_line_pnt_inside_polygon)
						{
							// create new point
							GEOLIB::Point* s (
							        bounding_polygon->
							        getIntersectionPointPolygonLine (
							                *(all_polylines[k])->
							                getPoint(
							                        j),
							                *(all_polylines[k])->
							                getPoint(j + 1)));
							if (s != NULL)
							{
								_out << "Point(" <<
								_n_pnt_offset << ") = {" <<
								(*s)[0] << ","
								<< (*s)[1] << "," <<
								(*s)[2] << "," <<
								mesh_density
								<<
								"}; // new end point of polyline "
								<< k << std::endl;
								// write line
								_out << "Line(" << _n_lines + j <<
								") = {" <<
								(all_polylines[k])->
								getPointID(j) << ","
								              << _n_pnt_offset <<
								"};" <<
								std::endl;
								// write line as constraint
								_out << "Line {" << _n_lines + j <<
								"} In Surface {" <<
								_n_plane_sfc - 1 << "};" <<
								std::endl;
								_n_pnt_offset++;
								delete s;
							}
						}
						if (!begin_line_pnt_inside_polygon &&
						    end_line_pnt_inside_polygon)
						{
							// create new point
							GEOLIB::Point* s (
							        bounding_polygon->
							        getIntersectionPointPolygonLine(
							                *(all_polylines[k])->
							                getPoint(
							                        j),
							                *(all_polylines[k])->
							                getPoint(j + 1)));
							if (s != NULL)
							{
								_out << "Point(" <<
								_n_pnt_offset << ") = {" <<
								(*s)[0] << ","
								<< (*s)[1] << "," <<
								(*s)[2] << "," <<
								mesh_density
								<<
								"}; // new end point of polyline "
								<< k << std::endl;
								// write line
								_out << "Line(" << _n_lines + j <<
								") = {" << _n_pnt_offset <<
								"," <<
								(all_polylines[k])->getPointID(
								        j + 1)
								<< "};" << std::endl;
								// write line as constraint
								_out << "Line {" << _n_lines + j <<
								"} In Surface {" <<
								_n_plane_sfc - 1 << "};" <<
								std::endl;
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

	// write stations as constraints
	_out << "// Stations" << std::endl;
	const size_t n_stations (all_stations.size());
	for (size_t k(0); k < n_stations; k++)
		if (bounding_polygon->isPntInPolygon (*(all_stations[k])))
		{
			_out << "Point(" << _n_pnt_offset + k << ") = {" <<
			(*(all_stations[k]))[0] << ","
			     << (*(all_stations[k]))[1] << "," << (*(all_stations[k]))[2] << "," <<
			mesh_density
			     << "};" << std::endl;
			_out << "Point {" << _n_pnt_offset + k << "} In Surface {" <<
			_n_plane_sfc - 1 << "};" << std::endl;
		}
	_n_pnt_offset += n_stations;

	std::cout << "ok" << std::endl;
}

void GMSHInterface::fetchGeometries (GEOLIB::GEOObjects const & geo,
                                     std::vector<std::string> const & selected_geometries,
                                     std::vector<GEOLIB::Point*>& all_points,
                                     std::vector<GEOLIB::Polyline*>& all_polylines,
                                     std::vector<GEOLIB::Point*>& all_stations) const
{
	// get names of all available data sources except stations
	std::vector<std::string> geo_names;
	// get station names
	std::vector<std::string> geo_station_names;

	for (std::vector<std::string>::const_iterator it (selected_geometries.begin());
	     it != selected_geometries.end(); ++it)
	{
		if ((geo.getPointVecObj (*it))->getType() == GEOLIB::PointVec::POINT)
			geo_names.push_back (*it);
		else if ((geo.getPointVecObj (*it))->getType() == GEOLIB::PointVec::STATION)
			geo_station_names.push_back (*it);
	}

	size_t pnt_offset (0);
	// fetch points and polylines and add them to the vectors, add points to the QuadTree
	for (std::vector<std::string>::const_iterator it (geo_names.begin());
	     it != geo_names.end(); ++it)
	{
		// get data from geo
#ifndef NDEBUG
		std::cout << "fetch geometrical data for " << *it << " " << std::flush;
#endif
		const std::vector<GEOLIB::Point*>* pnts (geo.getPointVec (*it));
		const std::vector<GEOLIB::Polyline*>* plys (geo.getPolylineVec (*it));
#ifndef NDEBUG
		std::cout << "ok" << std::endl;
#endif

		if (pnts)
		{
			// insert points into vector all_points
			all_points.insert (all_points.end(), pnts->begin(), pnts->end());

			if (plys)
				for (size_t k(0); k < plys->size(); k++)
				{
					size_t pos (all_polylines.size());
					// insert new polyline
					all_polylines.push_back (new GEOLIB::Polyline (all_points));
					// copy points
					for (size_t j(0); j < (*plys)[k]->getNumberOfPoints(); j++)
						// set points of polyline
						(all_polylines[pos])->addPoint (
						        pnt_offset +
						        ((*plys)[k])->
						        getPointID(j));
				}
			pnt_offset += pnts->size();
		}
	}

	for (std::vector<std::string>::const_iterator it (geo_station_names.begin());
	     it != geo_station_names.end(); ++it)
	{
		// get data from geo
#ifndef NDEBUG
		std::cout << "fetch station data for " << *it << " " << std::flush;
#endif
		const std::vector<GEOLIB::Point*>* pnts (geo.getPointVec (*it));
#ifndef NDEBUG
		std::cout << "ok" << std::endl;
#endif
		// insert points into vector all_stations
		all_stations.insert (all_stations.end(), pnts->begin(), pnts->end());
	}
}

GEOLIB::Polygon* GMSHInterface::getBoundingPolygon (
        std::vector<GEOLIB::Polyline*> const & all_polylines,
        size_t& bp_idx) const
{
	GEOLIB::Polygon* bounding_polygon(NULL);
	const size_t n_polylines (all_polylines.size());
	for (size_t k(0); k < n_polylines; k++)
		if (all_polylines[k]->isClosed ()) // == Polygon
		{
			if (bounding_polygon) // we have already a bounding polygon
			{
				if (!bounding_polygon->isPolylineInPolygon (*(all_polylines[k])))
				{
					GEOLIB::Polygon* tmp_polygon (new GEOLIB::Polygon (
					                                      *(all_polylines[k])));
					if (tmp_polygon->isPolylineInPolygon (*bounding_polygon))
					{
						// found new bounding polygon
						delete bounding_polygon;
						bounding_polygon = tmp_polygon;
						bp_idx = k;
					}
					else
						std::cerr <<
						"INFO: there is no inclusion relation between the polygons "
						          << k << " and " << bp_idx << std::endl;
				}
			}
			else
			{
				bounding_polygon = new GEOLIB::Polygon (*(all_polylines[k]));
				bp_idx = k;
			}
		}

	if (!bounding_polygon)
	{
		std::cerr <<
		"WARNING: GMSHInterface::writeAllDataToGMSHInputFile: did not found bounding polygon - abort writing"
		          << std::endl;
		return NULL;
	}

	return bounding_polygon;
}

void GMSHInterface::writeBoundingPolygon (GEOLIB::Polygon const* const bounding_polygon )
{
	std::cout << "write bounding polygon ... " << std::flush;
	// write bounding polygon
	size_t s (bounding_polygon->getNumberOfPoints());
	// write line segments (= Line) of the polyline
	for (size_t j(0); j < s - 1; j++)
		_out << "Line(" << _n_lines + j << ") = {" <<  _n_pnt_offset +
		bounding_polygon->getPointID(j) << ","
		     << _n_pnt_offset + bounding_polygon->getPointID(j + 1) << "};" << std::endl;
	// write the line segments contained in the polyline (=Line Loop)
	_out << "Line Loop (" << _n_lines + s - 1 << ") = {";
	for (size_t j(0); j < s - 2; j++)
		_out << _n_lines + j << ",";
	_out << _n_lines + s - 2 << "};" << std::endl;
	_n_lines += s;
	// write plane surface
	_out << "Plane Surface (" << _n_plane_sfc << ") = {" << _n_lines - 1 << "};" << std::endl;
	_n_plane_sfc++;
	std::cout << "ok" << std::endl;
}

void GMSHInterface::addPointsAsConstraints(const std::vector<GEOLIB::Point*> &points,
                                           const std::vector<GEOLIB::Polyline*> &polylines,
                                           std::map<size_t,size_t> geo2gmsh_surface_id_map,
                                           GEOLIB::QuadTree<GEOLIB::Point>* quad_tree)
{
	std::vector<GEOLIB::Polygon*> polygons;
	for (size_t j = 0; j < polylines.size(); j++)
		if (polylines[j]->isClosed())
		{
			GEOLIB::Polygon* pgn = new GEOLIB::Polygon(*polylines[j]);
			polygons.push_back(pgn);
			//geo2gmsh_surface_id_map[polygons.size()-1] = geo2gmsh_surface_id_map[j]; // this should be the same as above but do you wanna take the risk?
		}

	size_t nPoints = points.size();
	for (size_t i = 0; i < nPoints; i++)
	{
		std::list<size_t> surrounding_polygons;
		for (size_t j = 0; j < polygons.size(); j++)
			if (polygons[j]->isPntInPolygon(*(points[i])))
				surrounding_polygons.push_back(j);

		if (!surrounding_polygons.empty())
		{
			for (std::list<size_t>::iterator it = surrounding_polygons.begin();
			     it != surrounding_polygons.end(); ++it)
				for (std::list<size_t>::iterator jt = surrounding_polygons.begin();
				     jt != surrounding_polygons.end(); )
				{
					if (it != jt)
					{
						if (polygons[*it]->isPolylineInPolygon(*(polygons[*
						                                                  jt
						                                         ])))
							it = surrounding_polygons.erase(it);
						else
							++it;
					}
					else
						++jt;
				}

			_n_pnt_offset++;
			_out << "Point(" << _n_pnt_offset << ") = {" << (*points[i])[0] << "," <<
			(*points[i])[1] << "," << (*points[i])[2];
			if (quad_tree)
			{
				GEOLIB::Point ll, ur;
				quad_tree->getLeaf(*(points[i]), ll, ur);
				_out << "," << (0.3 * (ur[0] - ll[0]));
			}
			_out << "};" << std::endl;
			_out << "Point {" << _n_pnt_offset << "} In Surface {" <<
			geo2gmsh_surface_id_map[*(surrounding_polygons.begin())] << "};" <<
			std::endl;
		}
	}
}

std::vector<GEOLIB::Point*> GMSHInterface::getStationPoints(const GEOLIB::GEOObjects &geo)
{
	std::vector<GEOLIB::Point*> station_points;
	std::vector<std::string> stn_names;
	geo.getStationNames(stn_names);

	for (std::vector<std::string>::const_iterator it (stn_names.begin()); it != stn_names.end();
	     ++it)
	{
		const std::vector<GEOLIB::Point*>* pnts (geo.getPointVec (*it));
		station_points.insert (station_points.end(), pnts->begin(), pnts->end());
	}
	return station_points;
}

std::vector<GEOLIB::Point*> GMSHInterface::getSteinerPoints(
        GEOLIB::QuadTree<GEOLIB::Point>* quad_tree)
{
	std::vector<GEOLIB::Point*> steiner_points;
	std::list<GEOLIB::QuadTree<GEOLIB::Point>*> leaf_list;
	quad_tree->getLeafs (leaf_list);
	_out << "// Steiner points" << std::endl;
	for (std::list<GEOLIB::QuadTree<GEOLIB::Point>*>::const_iterator it (leaf_list.begin());
	     it != leaf_list.end(); it++)
		if ((*it)->getPoints().empty())
		{
			// compute point from square
			GEOLIB::Point ll, rr;
			(*it)->getSquarePoints (ll, rr);
			GEOLIB::Point* mid_point =
			        new GEOLIB::Point(0.5 * (rr[0] + ll[0]),
			                          0.5 * (rr[1] + ll[1]),
			                          0.5 * (rr[2] + ll[2]));
			steiner_points.push_back(mid_point);
		}
	return steiner_points;
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
					mesh->ele_vector[i]->nodes_index[j] =
					        gmsh_id[mesh->ele_vector[i]->GetNodeIndex(j) + 1];

			for (size_t i = 0; i < mesh->nod_vector.size(); i++)
				mesh->nod_vector[i]->SetIndex(i);
			// END OF: ordering nodes and closing gaps TK
		} /*End while*/
	}
	in.close();
}
} // end namespace FileIO
