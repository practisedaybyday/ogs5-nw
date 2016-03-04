/**
 * \file MSHEnums.h
 * 15/11/2010 KR initial implementation
 *
 */

#ifndef MSHENUMS_H
#define MSHENUMS_H

#include <string>

/**
 * \brief Types of mesh elements supported by OpenGeoSys.
 * Classification and associated int are identical with
 * identifier in MeshLib::CElem::geo_type
 */
struct MshElemType
{
	enum type {
		LINE = 1,
		QUAD = 2,
		QUAD8 = 3,
		HEXAHEDRON = 4,
		TRIANGLE = 5,
		TETRAHEDRON = 6,
		PRISM = 7,
		PYRAMID = 8,
		INVALID = -1,
		LAST = PYRAMID
	};
};

struct MshQualityType
{
	enum type {
		INVALID = 0,
		AREA,
		VOLUME,
		EDGERATIO,
		EQUIANGLESKEW
	};
};

/// Given a MshElemType this returns the appropriate string.
const std::string MshElemType2String(const MshElemType::type t);

/// Given a string describing an element type this returns the corresponding MshElemType.
MshElemType::type String2MshElemType(const std::string &s);

const std::string MshQualityType2String(const MshQualityType::type t);

#endif //MSHENUMS_H
