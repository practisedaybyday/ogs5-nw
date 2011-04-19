/**
 * \file FEMCondition.cpp
 * 25/11/2010 KR inital implementation
 *
 */

#include "FEMCondition.h"

#include "rf_bc_new.h"
#include "rf_ic_new.h"
#include "rf_st_new.h"

FEMCondition::FEMCondition(const std::string &geometry_name, CondType t) 
: _type(t), _geoObject(NULL), _geoName("[unspecified]"), _associated_geometry(geometry_name)
{
	this->setProcessType(INVALID_PROCESS);
	this->setProcessPrimaryVariable(INVALID_PV);
	this->setGeoType(GEOLIB::INVALID);
	this->setProcessDistributionType(FiniteElement::INVALID_DIS_TYPE);
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


BoundaryCondition::BoundaryCondition(const CBoundaryCondition &bc, const std::string &geometry_name)
: FEMCondition(geometry_name, FEMCondition::BOUNDARY_CONDITION)
{
	this->setProcessType(bc.getProcessType());
	this->setProcessPrimaryVariable(bc.getProcessPrimaryVariable());
	this->setGeoType(bc.getGeoType());
	this->setGeoName(bc.getGeoName());
	this->setProcessDistributionType(bc.getProcessDistributionType());

	if (this->getProcessDistributionType() == FiniteElement::CONSTANT) 
		this->setDisValue(bc.getGeoNodeValue());
	else if (this->getProcessDistributionType() == FiniteElement::LINEAR) 
		this->setLinearDisValues(bc.getPointsWithDistribedBC(), bc.getDistribedBC());
}

InitialCondition::InitialCondition(const CInitialCondition &ic, const std::string &geometry_name)
: FEMCondition(geometry_name, FEMCondition::INITIAL_CONDITION)
{
	this->setProcessType(ic.getProcessType());
	this->setProcessPrimaryVariable(ic.getProcessPrimaryVariable());
	this->setGeoType(ic.getGeoType());
	std::string geo_name = (ic.getGeoType() == GEOLIB::GEODOMAIN) ? "Domain" : ic.getGeoName();
	this->setGeoName(geo_name);
	this->setProcessDistributionType(ic.getProcessDistributionType());

	if (this->getProcessDistributionType() == FiniteElement::CONSTANT)
		this->setDisValue(ic.getGeoNodeValue());
}

SourceTerm::SourceTerm(const CSourceTerm &st, const std::string &geometry_name)
: FEMCondition(geometry_name, FEMCondition::SOURCE_TERM)
{
	this->setProcessType(st.getProcessType());
	this->setProcessPrimaryVariable(st.getProcessPrimaryVariable());
	this->setGeoType(st.getGeoType());
	this->setGeoName(st.getGeoName());
	this->setProcessDistributionType(st.getProcessDistributionType());

	if (this->getProcessDistributionType() == FiniteElement::CONSTANT)
		this->setDisValue(st.getGeoNodeValue());
	else if (this->getProcessDistributionType() == FiniteElement::LINEAR) 
		this->setLinearDisValues(st.getPointsWithDistribedST(), st.getDistribedST());

}
