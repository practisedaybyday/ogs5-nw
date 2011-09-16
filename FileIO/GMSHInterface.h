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

namespace FileIO {

/**
 * \brief Reads and writes GMSH-files to and from OGS data structures.
 */
class GMSHInterface {
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
	 * writes the geometric data (Points, Polylines, Surfaces) into a file of the GMSH format
	 * @param proj_name Name of the geometry that will be included in the geo-file
	 * @param geo Container for geometric information
	 * @parem useStationsAsContraints If true, station data is included as constraints for meshing of surfaces (via addStationsAsConstraints()).
	 * @return if the file stream can be opened the method returns true, else it returns false
	 */
	bool writeGMSHInputFile(const std::string &proj_name, const GEOLIB::GEOObjects& geo, bool useStationsAsContraints = false);

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
			double mesh_density_scaling = 0.3, double mesh_density_scaling_station_pnts = 0.05);

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

	void writeGMSHPoints(const std::vector<GEOLIB::Point*>& pnt_vec);
	void writeGMSHPolyline (const GEOLIB::Polyline* ply, const size_t offset);
	void writeGMSHPolylines(const std::vector<GEOLIB::Polyline*>& ply_vec);
	size_t writeGMSHPolygon(const GEOLIB::Polygon& polygon, const size_t offset);
	void writePlaneSurface (std::list<size_t> const & polygon_list);

	static bool isGMSHMeshFile (const std::string& fname);

private:
	void fetchGeometries (GEOLIB::GEOObjects const & geo,
			std::vector<std::string> const & selected_geometries,
			std::vector<GEOLIB::Point*>& all_points,
			std::vector<GEOLIB::Polyline*>& all_polylines,
			std::vector<GEOLIB::Point*>& all_stations) const;
	std::list<size_t> findHolesInsidePolygon(const std::vector<GEOLIB::Polyline*> *plys, size_t i, std::map<size_t,size_t> geo2gmsh_polygon_id_map);
	GEOLIB::Polygon* getBoundingPolygon (std::vector<GEOLIB::Polyline*> const & all_polylines, size_t &bp_idx) const;
	void writeBoundingPolygon (GEOLIB::Polygon const * const bounding_polygon );
	/// Add station points as constraints to the current geometry
	void addStationsAsConstraints(const std::string &proj_name, const GEOLIB::GEOObjects& geo, std::map<size_t,size_t> geo2gmsh_surface_id_map);
	size_t _n_pnt_offset;
	size_t _n_lines;
	size_t _n_plane_sfc;
	std::ofstream _out;
	//std::list<size_t> _polygon_list;
};

}

#endif /* GMSHINTERFACE_H_ */
