/**
 * \file ConversionTools.h
 * 2012/04/11 KR Initial implementation
 */

#ifndef CONVERSIONTOOLS_H
#define CONVERSIONTOOLS_H

#include <vector>
#include <QString>

#include "FEMCondition.h"

class ConversionTools
{
public:
	static void getFEMConditionsFromASCIIFile(const QString &file_name, std::vector<FEMCondition*> &conditions);

private:
	static std::ios::pos_type readBoundaryCondition(std::vector<FEMCondition*> &conditions, std::ifstream &in, GEOLIB::GEOObjects &geo_objects, std::string geo_name);
	static std::ios::pos_type readInitialCondition (std::vector<FEMCondition*> &conditions, std::ifstream &in, GEOLIB::GEOObjects &geo_objects, std::string geo_name);
	static std::ios::pos_type readSourceTerm       (std::vector<FEMCondition*> &conditions, std::ifstream &in, GEOLIB::GEOObjects &geo_objects, std::string geo_name);

};

#endif //CONVERSIONTOOLS_H
