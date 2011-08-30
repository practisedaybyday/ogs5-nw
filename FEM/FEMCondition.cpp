/**
 * \file FEMCondition.cpp
 * 25/11/2010 KR inital implementation
 *
 */

#include "FEMCondition.h"
#include "ProcessInfo.h"

#include "GEOObjects.h" //for SourceTerm
#include "GridAdapter.h"

FEMCondition::FEMCondition(const std::string &geometry_name, CondType t) 
: _type(t), _geoObject(NULL), _geoName("[unspecified]"), _associated_geometry(geometry_name)
{
	this->setProcessType(INVALID_PROCESS);
	this->setProcessPrimaryVariable(INVALID_PV);
	this->setGeoType(GEOLIB::INVALID);
	this->setProcessDistributionType(FiniteElement::INVALID_DIS_TYPE);
}

FEMCondition::FEMCondition(const std::string &geometry_name, ProcessType pt, PrimaryVariable pv, GEOLIB::GEOTYPE gt, const std::string &gn, FiniteElement::DistributionType dt, CondType ct)
: ProcessInfo(pt, pv, NULL), GeoInfo(gt, NULL), DistributionInfo(dt), _type(ct), _geoObject(NULL), _geoName(gn), _associated_geometry(geometry_name)
{
}

std::string FEMCondition::condTypeToString(CondType type)
{
	if (type==FEMCondition::BOUNDARY_CONDITION) return "Boundary Conditions";
	else if (type==FEMCondition::INITIAL_CONDITION) return "Initial Conditions";
	else if (type==FEMCondition::SOURCE_TERM) return "Source Terms";
	else return "Unspecified";
}

void FEMCondition::setLinearDisValues(const std::vector<int> &point_ids, const std::vector<double> &point_values)
{
	for (size_t i=0; i<point_ids.size(); i++)
	{
		this->_disValue.push_back(point_ids[i]);
		this->_disValue.push_back(point_values[i]);
	}
}

