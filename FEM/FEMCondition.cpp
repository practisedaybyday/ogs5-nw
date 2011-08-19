/**
 * \file FEMCondition.cpp
 * 25/11/2010 KR inital implementation
 *
 */

#include "FEMCondition.h"

#include "rf_bc_new.h"
#include "rf_ic_new.h"
#include "rf_st_new.h"
#include "ProcessInfo.h"

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
/*
	else if (this->getProcessDistributionType() == FiniteElement::DIRECT) 
	{
		this->setGeoType(GEOLIB::POINT);
		this->setGeoName("Point");
		this->setDisValue(st.getGeoNodeValue());

		st.getProcess()->st_node_value[0]->geo_node_number;
		this->setGeoObj(
		st.getProcess()->st_node_value[0]->node_value;
		//st.DirectAssign();
	}
*/
	else std::cout << "Error in SourceTerm() - Unknown Process Distribution Type \"" << FiniteElement::convertDisTypeToString(st.getProcessDistributionType()) << "\"..." << std::endl;
}

std::vector<FEMCondition*> SourceTerm::createDirectSourceTerms(const std::vector<CSourceTerm*> &st_vector, const std::string &geo_name, const std::vector<GEOLIB::Point*> *new_points)
{
	// read source term file and make sure it's really DIRECT-STs
	std::vector<FEMCondition*> conditions;
	
	for (std::vector<CSourceTerm*>::const_iterator it = st_vector.begin(); it != st_vector.end(); ++it)
	{
		if ((*it)->getProcessDistributionType() == FiniteElement::DIRECT)
		{
			std::vector< std::pair<size_t, double> > node_values;
			SourceTerm::getDirectNodeValues((*it)->fname, node_values);

			// create one constant boundary condition for every point specified in the DIRECT-array
			size_t nNodes = node_values.size();
			for (size_t i=0; i<nNodes; i++)
			{
				std::stringstream out;
				out << node_values[i].first;
				SourceTerm* st = new SourceTerm(geo_name);
				st->setProcessType((*it)->getProcessType());
				st->setProcessPrimaryVariable((*it)->getProcessPrimaryVariable());
				st->setGeoType(GEOLIB::POINT);
				st->setGeoName(out.str());
				st->setProcessDistributionType(FiniteElement::CONSTANT);
				st->setDisValue( node_values[i].second );
				conditions.push_back(st);
			}
		}
		else 
			std::cout << "Error: no DIRECT distribution type" << std::endl;
	}
	return conditions;
}

void SourceTerm::getDirectNodeValues(const std::string &filename, std::vector< std::pair<size_t, double> > &node_values)
{
	std::ifstream in(filename.c_str());
	if (!in.is_open())
	{
		std::cout << "Error in getNodeValues() - Could not find file for direct node values..." << std::endl;
		return;
	}

	std::stringstream str_in;
	std::string line("");
	size_t idx(0);
	double val(0);

	while ( getline(in, line) )
	{
		if (line.find("#STOP") != std::string::npos) return;
		str_in << line;
		str_in >> idx >> val;
		node_values.push_back(std::pair<size_t, double>(idx, val));
		str_in.clear();
	}
}