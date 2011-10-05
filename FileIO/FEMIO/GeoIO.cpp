/*
 * GeoIO.cpp
 *
 *  Created on: Sep 29, 2010
 *      Author: TF
 */

#include <sstream>

#include "FEMIO/GeoIO.h"

// FEM
#include "readNonBlankLineFromInputStream.h"

namespace FileIO {

bool GeoIO::readGeoInfo (GeoInfo* geo_info, std::istream& in_str, std::string& geo_name,
			const GEOLIB::GEOObjects& geo_obj, const std::string& unique_geo_name)
{
	std::stringstream strstream;
	strstream.str(readNonBlankLineFromInputStream(in_str));
	std::string geo_type_name;
	strstream >> geo_type_name;

	if (geo_type_name.find("POINT") != std::string::npos) {
		strstream >> geo_name;
		const GEOLIB::Point *pnt(
				(geo_obj.getPointVecObj(unique_geo_name))->getElementByName(geo_name));
		if (pnt == NULL) {
			std::cerr << "ERROR in GeoIO::readGeoInfo: point name \"" << geo_name
					<< "\" not found!" << std::endl;
#ifdef OGS_USE_QT
			return false;
#else
			exit(1);
#endif
		}
		geo_info->setGeoType(GEOLIB::POINT);
		geo_info->setGeoObj(pnt);
		strstream.clear();
		return true;
	}

	if (geo_type_name.find("POLYLINE") != std::string::npos) {
		geo_info->setGeoType(GEOLIB::POLYLINE);
		strstream >> geo_name;
		const GEOLIB::Polyline *ply(
				(geo_obj.getPolylineVecObj(unique_geo_name))->getElementByName(geo_name));
		if (ply == NULL) {
			std::cerr << "error in GeoIO::readGeoInfo: polyline name \"" << geo_name
					<< "\" not found!" << std::endl;
#ifdef OGS_USE_QT
			return false;
#else
			exit(1);
#endif
		}
		geo_info->setGeoObj(ply);
		strstream.clear();
		return true;
	}

	if (geo_type_name.find("SURFACE") != std::string::npos) {
		geo_info->setGeoType(GEOLIB::SURFACE);
		strstream >> geo_name;
		GEOLIB::SurfaceVec const* sfc_vec (geo_obj.getSurfaceVecObj(unique_geo_name));
		if (sfc_vec) {
			const GEOLIB::Surface *sfc(sfc_vec->getElementByName(geo_name));
			if (sfc == NULL) {
				std::cerr << "error in GeoIO::readGeoInfo: surface name \"" << geo_name
					<< "\" not found!" << std::endl;
#ifdef OGS_USE_QT
				return false;
#else
				exit(1);
#endif
			}
			geo_info->setGeoObj(sfc);
		} else {
			std::cerr << "error in GeoIO::readGeoInfo: surface vector not found!" << std::endl;
#ifdef OGS_USE_QT
			return false;
#else
			exit(1);
#endif
		}
		strstream.clear();
		return true;
	}

	if (geo_type_name.find("VOLUME") != std::string::npos) {
		geo_info->setGeoType(GEOLIB::VOLUME);
		strstream >> geo_name;
		strstream.clear();
		return true;
	}

	if (geo_type_name.find("DOMAIN") != std::string::npos) {
		geo_info->setGeoType(GEOLIB::GEODOMAIN);
		strstream >> geo_name;
		strstream.clear();
		return true;
	}

	return false;
}

}
