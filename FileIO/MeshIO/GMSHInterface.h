/*
 * GMSHInterface.h
 *
 *  Created on: Apr 29, 2010
 *      Author: TF
 */

#ifndef GMSHINTERFACE_H_
#define GMSHINTERFACE_H_

#include <string>

// GEOLIB
#include "GEOObjects.h"
#include "Polygon.h"
//#include "Point.h"
#include "QuadTree.h"

namespace MeshLib
{
class CFEMesh;
}

namespace FileIO
{
/**
 * \brief Reads and writes GMSH-files to and from OGS data structures.
 */
class GMSHInterface
{
public:
	/**
	 * constructor opens a file stream with the given file name
	 * @param fname file name
	 * @return
	 */
	GMSHInterface (const std::string &fname);
	/**
	 * destructor closes the stream
	 * @return
	 */
	~GMSHInterface ();
	/**
	 * writes the geometric data (Points, Polylines, Surfaces) into a file of the GMSH geometric format
	 * @param proj_name Name of the geometry that will be included in the geo-file
	 * @param geo Container for geometric information
	 * @param useStationsAsContraints If true, station data is included as constraints for meshing of surfaces (via addStationsAsConstraints()).
	 * @param useSteinerPoints If true, additional points will be generated based on a quadtree such that a certained pre-defined point-density is set.
	 * @return if the file stream can be opened the method returns true, else it returns false
	 */
	bool writeGMSHInputFile(const std::string &proj_name,
	                        const GEOLIB::GEOObjects& geo,
	                        bool useStationsAsContraints = false,
							bool useSteinerPoints = false);

	/**
	 * Method writes selected data to the stream (opened from constructor) employing a Quadtree for
	 * adaptive mesh generation.
	 *
	 * @param geo object managing all geometric information
	 * @param selected_geometries geometric information that should be included into the mesh process
	 * @param number_of_point_per_quadtree_node maximum number of points per Quadtree leaf
	 * (see class \sa Quadtree for documentation)
	 * @param mesh_density_scaling The mesh density at a point depends on the edge size
	 * of the Quadtree leaf the point is located in. The mesh_density is the edge size
	 * multiplied with the scaling factor mesh_density_scaling.
	 * @param mesh_density_scaling_station_pnts The mesh density at a station depends on the edge size
	 * of the Quadtree leaf the point is located in. The mesh_density is the edge size
	 * multiplied with the scaling factor mesh_density_scaling_station_pnts.
	 */
	void writeAllDataToGMSHInputFile (GEOLIB::GEOObjects& geo,
	                                  std::vector<std::string> const & selected_geometries,
	                                  size_t number_of_point_per_quadtree_node = 10,
	                                  double mesh_density_scaling = 0.3,
	                                  double mesh_density_scaling_station_pnts = 0.05);

	/**
	 * Method writes selected data to the stream (opened from constructor) for mesh generation.
	 *
	 * @param geo object managing all geometric information
	 * @param selected_geometries geometric information that should be included into the mesh process
	 * @param mesh_density The mesh density at a point.
	 */
	void writeAllDataToGMSHInputFile (GEOLIB::GEOObjects& geo,
	                                  std::vector<std::string> const & selected_geometries,
	                                  double mesh_density);

	std::string writeAllDataToGMSHInput (GEOLIB::GEOObjects& geo,
	                                    std::vector<std::string> const & selected_geometries,
	                                    size_t number_of_point_per_quadtree_node = 10,
	                                    double mesh_density_scaling = 0.3,
	                                    double mesh_density_scaling_station_pnts = 0.05);

	void writeGMSHPoints(const std::vector<GEOLIB::Point*>& pnt_vec,
	                     GEOLIB::QuadTree<GEOLIB::Point>* quad_tree = NULL);
	void writeGMSHPolyline (std::vector<GEOLIB::Polyline*> const& plys, size_t ply_id, size_t offset);
	void writeGMSHPolylines(const std::vector<GEOLIB::Polyline*>& ply_vec);
	size_t writeGMSHPolygon(std::vector<GEOLIB::Polyline*> const& ply_vec, size_t ply_id, size_t offset);
	void writePlaneSurface (std::list<size_t> const & polygon_list, bool respect_holes);

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

private:
	GEOLIB::QuadTree<GEOLIB::Point>* createQuadTreeFromPoints(
	        std::vector<GEOLIB::Point*> points,
	        size_t number_of_point_per_quadtree_node);
	void fetchGeometries (GEOLIB::GEOObjects const & geo,
	                      std::vector<std::string> const & selected_geometries,
	                      std::vector<GEOLIB::Point*>& all_points,
	                      std::vector<GEOLIB::Polyline*>& all_polylines,
	                      std::vector<GEOLIB::Point*>& all_stations) const;
	std::list<size_t> findHolesInsidePolygon(const std::vector<GEOLIB::Polyline*>* plys, size_t i) const;
	/**
	 * Method tests if the i-th polygon is included in any other polygon.
	 * @param plys the vector of polylines (closed polylines are considered as polygon)
	 * @param i the polygon for the test
	 * @param j if the method returns true, the j-th polygon contains the i-th polygon
	 * @return true, if polygon is in an other polygon included, else false
	 */
	bool isPolygonInOtherPolygon(const std::vector<GEOLIB::Polyline*>* plys, size_t i, size_t &j) const;
	GEOLIB::Polygon* getBoundingPolygon (std::vector<GEOLIB::Polyline*> const & all_polylines,
	                                     size_t &bp_idx) const;
	std::vector<GEOLIB::Point*> getStationPoints(const GEOLIB::GEOObjects &geo_objects);
	std::vector<GEOLIB::Point*> getSteinerPoints(GEOLIB::QuadTree<GEOLIB::Point>* quad_tree);
	void writeBoundingPolygon (GEOLIB::Polygon const* const bounding_polygon );
	/// Adds a point-array (such as stations) as constraints to the geometry given by proj_name.
	void addPointsAsConstraints(const std::vector<GEOLIB::Point*> &points,
	                            const std::vector<GEOLIB::Polyline*> &polylines,
	                            std::map<size_t,size_t> geo2gmsh_surface_id_map,
	                            GEOLIB::QuadTree<GEOLIB::Point>* quad_tree = NULL,
	                            double mesh_density = 0.5);

	size_t _n_pnt_offset;
	size_t _n_lines;
	size_t _n_plane_sfc;
	std::stringstream _out;
	std::ofstream _out_file;

	/**
	 * this private class is for storing meta data, i.e. the mapping between
	 * polyline id, the number of the first and the last line segments
	 * in gmsh file that are describing the polyline in the gmsh file
	 */
	struct PolylineGMSHMapping {
		PolylineGMSHMapping (size_t ply_id, size_t start_id, size_t end_id) :
			_ply_id (ply_id), _gmsh_line_start_id(start_id), _gmsh_line_end_id(end_id)
		{}
		size_t _ply_id;
		size_t _gmsh_line_start_id;
		size_t _gmsh_line_end_id;
	};
	/**
	 * vector contains the mapping between polyline ids and first and last segment in gmsh file
	 */
	std::vector<PolylineGMSHMapping> _ply_id_gmsh_line_mapping;
};
}

#endif /* GMSHINTERFACE_H_ */
