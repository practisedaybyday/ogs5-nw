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


BoundaryCondition::BoundaryCondition(const CBoundaryCondition &bc, const std::string &geometry_name)
: FEMCondition(geometry_name, bc.getProcessType(), bc.getProcessPrimaryVariable(), bc.getGeoType(), bc.getGeoName(), 
			   bc.getProcessDistributionType(), FEMCondition::BOUNDARY_CONDITION)
{
	if (this->getProcessDistributionType() == FiniteElement::CONSTANT || this->getProcessDistributionType() == FiniteElement::CONSTANT_NEUMANN) 
		this->setDisValue(bc.getGeoNodeValue());
	else if (this->getProcessDistributionType() == FiniteElement::LINEAR || this->getProcessDistributionType() == FiniteElement::LINEAR_NEUMANN) 
		this->setLinearDisValues(bc.getPointsWithDistribedBC(), bc.getDistribedBC());
}

InitialCondition::InitialCondition(const CInitialCondition &ic, const std::string &geometry_name)
: FEMCondition(geometry_name, ic.getProcessType(), ic.getProcessPrimaryVariable(), ic.getGeoType(), (ic.getGeoType() == GEOLIB::GEODOMAIN) ? "Domain" : ic.getGeoName(), 
			   ic.getProcessDistributionType(), FEMCondition::INITIAL_CONDITION)
{
	if (this->getProcessDistributionType() == FiniteElement::CONSTANT)
		this->setDisValue(ic.getGeoNodeValue());
}

SourceTerm::SourceTerm(const CSourceTerm &st, const std::string &geometry_name)
: FEMCondition(geometry_name, st.getProcessType(), st.getProcessPrimaryVariable(), st.getGeoType(), st.getGeoName(), 
			   st.getProcessDistributionType(), FEMCondition::SOURCE_TERM)
{
	if (this->getProcessDistributionType() == FiniteElement::CONSTANT || this->getProcessDistributionType() == FiniteElement::CONSTANT_NEUMANN) 
		this->setDisValue(st.getGeoNodeValue());
	else if (this->getProcessDistributionType() == FiniteElement::LINEAR || this->getProcessDistributionType() == FiniteElement::LINEAR_NEUMANN) 
		this->setLinearDisValues(st.getPointsWithDistribedST(), st.getDistribedST());
}
