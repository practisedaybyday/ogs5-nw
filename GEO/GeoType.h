/*
 * GeoType.h
 *
 *  Created on: Jun 17, 2010
 *      Author: TF
 */

#ifndef GEOTYPE_H_
#define GEOTYPE_H_



namespace GEOLIB {

/**
 * \ingroup GEOLIB
 */

enum GEOTYPE {
	POINT,   //!< POINT
	POLYLINE,//!< POLYLINE
	SURFACE, //!< SURFACE
	VOLUME,  //!< VOLUME
	GEODOMAIN//!< GEODOMAIN
};

} // end namespace GEOLIB

#endif /* GEOTYPE_H_ */
