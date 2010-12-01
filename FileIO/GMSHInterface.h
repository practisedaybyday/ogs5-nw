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
	 * @param proj_name
	 * @param geo
	 * @return if the file stream can be opened the method returns true, else it returns false
	 */
	bool writeGMSHInputFile(const std::string &proj_name, const GEOLIB::GEOObjects& geo);

	void writeAllDataToGMSHInputFile (GEOLIB::GEOObjects& geo);

	void writeGMSHPoints(const std::vector<GEOLIB::Point*>& pnt_vec);
	void writeGMSHPolyline (const GEOLIB::Polyline* ply);
	void writeGMSHPolylines(const std::vector<GEOLIB::Polyline*>& ply_vec);
	void writeGMSHPolygon(const GEOLIB::Polygon& polygon);
	void writePlaneSurface ();

	static bool isGMSHMeshFile (const std::string& fname);

private:
	size_t _n_pnt_offset;
	size_t _n_lines;
	size_t _n_plane_sfc;
	std::ofstream _out;
	std::list<size_t> _polygon_list;
};

}

#endif /* GMSHINTERFACE_H_ */
