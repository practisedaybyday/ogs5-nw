/*
 * GMSHInterface.h
 *
 *  Created on: Apr 29, 2010
 *      Author: TF
 */

#ifndef GMSHINTERFACE_H_
#define GMSHINTERFACE_H_

#include <string>
#include <list>

// FileIO
#include "Writer.h"
#include "GMSHPoint.h"
#include "GMSHLineLoop.h"
#include "GMSHPolygonTree.h"
#include "GMSHMeshDensityStrategy.h"

// GEOLIB
#include "GEOObjects.h"
#include "Polygon.h"

namespace MeshLib
{
class CFEMesh;
}

namespace FileIO
{

namespace GMSH {

enum MeshDensityAlgorithm {
	NoMeshDensity = 0, //!< do not set the parameter
	FixedMeshDensity, //!< set the parameter with a fixed value
	AdaptiveMeshDensity //!< computing the mesh density employing a QuadTree
};

}

/**
 * \brief Reads and writes GMSH-files to and from OGS data structures.
 */
class GMSHInterface : public Writer
{
public:

	/**
	 *
	 * @param geo_objs reference tp instance of class GEOObject that maintains the geometries.
	 * 	The instance is used for preparation geometries for writing them to the gmsh file format.
	 * @param include_stations_as_constraints switch to enable writing stations as constraints
	 * @param mesh_density_algorithm one of the mesh density algorithms (\@see enum MeshDensityAlgorithm)
	 * @param param1 parameter that can be used for the mesh density algorithm
	 * @param param2 parameter that can be used for the mesh density algorithm
	 * @param param3 parameter that can be used for the mesh density algorithm
	 * @param selected_geometries vector of names of geometries, that should be employed for mesh generation.
	 * @return
	 */
	GMSHInterface (GEOLIB::GEOObjects & geo_objs,
					bool include_stations_as_constraints,
					GMSH::MeshDensityAlgorithm mesh_density_algorithm,
					double param1, double param2, size_t param3,
					std::vector<std::string> & selected_geometries);

	/**
	 * checks if there is a GMSH mesh file header
	 * @param fname the file name of the mesh (including the path)
	 * @return true, if the file seems to be a valid GMSH file, else false
	 */
	static bool isGMSHMeshFile (const std::string& fname);
	/**
	 * reads a mesh created by GMSH - this implementation is based on the former function GMSH2MSH
	 * @param fname the file name of the mesh (including the path)
	 * @param mesh the new mesh
	 * @return
	 */
	static void readGMSHMesh (std::string const& fname, MeshLib::CFEMesh* mesh);

protected:
	int write(std::ostream& stream);

private:
	/**
	 * 1. get and merge data from _geo_objs
	 * 2. compute topological hierarchy
	 * @param out
	 */
	void writeGMSHInputFile(std::ostream & out);

	void writePoints(std::ostream& out) const;
//	void writeLines(std::ostream& out) const;

	/**
	 * create vector of points
	 * @param bounding_polygon
	 * @param pnt_vec
	 * @param station
	 * @return
	 */
//	size_t createGMSHPointsInsideBoundingPolygon(GEOLIB::Polygon const*const bounding_polygon, const std::vector<GEOLIB::Point*> &pnt_vec, bool station=false);


//	std::list<size_t> findHolesInsidePolygon(const std::vector<GEOLIB::Polyline*>* plys, size_t i) const;
	/**
	 * Method tests if the i-th polygon is included in any other polygon.
	 * @param plys the vector of polylines (closed polylines are considered as polygon)
	 * @param i the polygon for the test
	 * @param j if the method returns true, the j-th polygon contains the i-th polygon
	 * @return true, if polygon is in an other polygon included, else false
	 */
//	bool isPolygonInOtherPolygon(const std::vector<GEOLIB::Polyline*>* plys, size_t i, size_t &j) const;

//	GEOLIB::Polygon* getBoundingPolygon (std::vector<GEOLIB::Polyline*> const & all_polylines,
//	                                     size_t &bp_idx) const;
//	void getStationPoints(const GEOLIB::GEOObjects &geo_objects, std::vector<GEOLIB::Point*>& station_points) const;
//	void getSteinerPoints(GEOLIB::QuadTree<GEOLIB::Point>* quad_tree, std::vector<GEOLIB::Point*>& steiner_pnts) const;
//	void createBoundingPolygon(GEOLIB::Polygon const* const bounding_polygon, size_t pnt_offset);

	/// Adds a point-array (such as stations) as constraints to the geometry given by proj_name.
//	void addPointsAsConstraints(const std::vector<GEOLIB::Point*> &points,
//					const std::vector<GEOLIB::Polyline*> &polylines, std::ostream &out);

//	void createPolylinesAsConstraints(GEOLIB::Polygon * bounding_polygon, size_t bp_idx,
//					std::vector<GEOLIB::Polyline*> const& all_polylines, std::ostream & out);

	size_t _n_lines;
	size_t _n_plane_sfc;

	GEOLIB::GEOObjects & _geo_objs;
	std::vector<std::string>& _selected_geometries;
	std::string _gmsh_geo_name;
	std::string _gmsh_station_name;
	std::list<GMSHPolygonTree*> _polygon_tree_list;

	bool _include_stations_as_constraints;

//	std::vector<size_t> _gmsh_pnt_id_to_geo_pnt_id_map;
	std::vector<FileIO::GMSHPoint*> _gmsh_pnts;
	std::vector<FileIO::GMSHLineLoop*> _gmsh_line_loops;

	GMSHMeshDensityStrategy *_mesh_density_strategy;

//	/**
//	 * this private class is for storing meta data, i.e. the mapping between
//	 * polyline id, the number of the first and the last line segments
//	 * in gmsh file that are describing the polyline in the gmsh file
//	 */
//	struct PolylineGMSHMapping {
//		PolylineGMSHMapping (size_t ply_id, size_t start_id, size_t end_id) :
//			_ply_id (ply_id), _gmsh_line_start_id(start_id), _gmsh_line_end_id(end_id)
//		{}
//		size_t _ply_id;
//		size_t _gmsh_line_start_id;
//		size_t _gmsh_line_end_id;
//	};
//	/**
//	 * vector contains the mapping between polyline ids and first and last segment in gmsh file
//	 */
//	std::vector<PolylineGMSHMapping> _ply_id_gmsh_line_mapping;
};
}

#endif /* GMSHINTERFACE_H_ */
