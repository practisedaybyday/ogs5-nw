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
	: _type(t), _geoName("[unspecified]"), _associated_geometry(geometry_name)
{
	this->setProcessType(FiniteElement::INVALID_PROCESS);
	this->setProcessPrimaryVariable(FiniteElement::INVALID_PV);
	this->setGeoType(GEOLIB::INVALID);
	this->setProcessDistributionType(FiniteElement::INVALID_DIS_TYPE);
}

FEMCondition::FEMCondition(const std::string &geometry_name, FiniteElement::ProcessType pt,
				FiniteElement::PrimaryVariable pv, GEOLIB::GEOTYPE gt, const std::string &gn,
				FiniteElement::DistributionType dt, CondType ct)
	: ProcessInfo(pt, pv, NULL),  GeoInfo(gt, NULL), DistributionInfo(dt), _type(ct),
	  _geoName(gn), _associated_geometry(geometry_name)
{
}

FEMCondition::FEMCondition(const FEMCondition &cond, CondType t)
	: ProcessInfo(cond.getProcessType(), cond.getProcessPrimaryVariable(), NULL),
	  GeoInfo(cond.getGeoType(), cond.getGeoObj()),
	  DistributionInfo(cond.getProcessDistributionType()),
	  _type(t),
	  _geoName(cond.getGeoName()),
	  _disValue(cond.getDisValue()),
	  _associated_geometry(cond.getAssociatedGeometryName())
{
}

std::string FEMCondition::condTypeToString(CondType type)
{
	if (type == FEMCondition::BOUNDARY_CONDITION)
		return "Boundary Conditions";
	else if (type == FEMCondition::INITIAL_CONDITION)
		return "Initial Conditions";
	else if (type == FEMCondition::SOURCE_TERM)
		return "Source Terms";
	else
		return "Unspecified";
}

void FEMCondition::setLinearDisValues(const std::vector< std::pair<size_t, double> > &dis_values)
{
	for (size_t i = 0; i < dis_values.size(); i++)
	{
		this->_disValue.push_back(dis_values[i].first);
		this->_disValue.push_back(dis_values[i].second);
	}
}

std::vector< std::pair<size_t, double> > FEMCondition::getDistributedPairs(std::vector<int> point_ids, std::vector<double> point_values)
{
	if (point_ids.size() == point_values.size())
	{
		size_t nValues (point_ids.size());
		std::vector< std::pair<size_t, double> > dis_values(nValues);
		for (size_t i=0; i<nValues; i++)
			dis_values.push_back( std::pair<size_t, double>(point_ids[i],point_values[i]) );
		return dis_values;
	}
	std::cout << "Error in SourceTerm() - size of linear distribution arrays doesn't match..." << std::endl;
	std::vector< std::pair<size_t, double> > dis_values(0);
	return dis_values;
}
