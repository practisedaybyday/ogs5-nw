/**
 * \file SourceTerm.cpp
 * 2011/08/30 KR inital implementation
 *
 */

#include "SourceTerm.h"
#include "rf_st_new.h"

SourceTerm::SourceTerm(const CSourceTerm &st, const std::string &geometry_name)
	: FEMCondition(geometry_name, st.getProcessType(), st.getProcessPrimaryVariable(),
	               st.getGeoType(), st.getGeoName(),
	               st.getProcessDistributionType(), FEMCondition::SOURCE_TERM)
{
	if (this->getProcessDistributionType() == FiniteElement::CONSTANT ||
	    this->getProcessDistributionType() == FiniteElement::CONSTANT_NEUMANN)
		this->setDisValue(st.getGeoNodeValue());
	else if (this->getProcessDistributionType() == FiniteElement::LINEAR ||
	         this->getProcessDistributionType() == FiniteElement::LINEAR_NEUMANN)
		this->setDisValues(this->getDistributedPairs(st.getPointsWithDistribedST(), st.getDistribedST()));
	else if (this->getProcessDistributionType() == FiniteElement::LINEAR)
	{
		this->_direct_file_name = st.fname;
	}
	else
		std::cout << "Error in SourceTerm() - Unknown Process Distribution Type \"" <<
		FiniteElement::convertDisTypeToString(st.getProcessDistributionType()) <<
		"\"..." <<
		std::endl;
}

std::vector<FEMCondition*> SourceTerm::createDirectSourceTerms(
        const std::vector<CSourceTerm*> &st_vector,
        const std::string &geo_name,
		const std::string &file_path)
{
	// read source term file and make sure it's really DIRECT-STs
	std::vector<FEMCondition*> conditions;

	for (std::vector<CSourceTerm*>::const_iterator it = st_vector.begin(); it != st_vector.end();
	     ++it)
	{
		if ((*it)->getProcessDistributionType() == FiniteElement::DIRECT)
		{
			std::vector< std::pair<size_t, double> > node_values;
			SourceTerm::getDirectNodeValues(file_path + "/" + (*it)->fname, node_values);

			SourceTerm* st = new SourceTerm(geo_name);
			st->setProcessType((*it)->getProcessType());
			st->setProcessPrimaryVariable((*it)->getProcessPrimaryVariable());
			st->setGeoType(GEOLIB::INVALID);
			st->setGeoName((*it)->fname);
			st->setProcessDistributionType(FiniteElement::DIRECT);
			st->setDisValues(node_values);
			conditions.push_back(st);
		}
		else
			std::cout << "Error: no DIRECT distribution type" << std::endl;
	}
	return conditions;
}

void SourceTerm::getDirectNodeValues(const std::string &filename,
                                     std::vector< std::pair<size_t, double> > &node_values)
{
	std::ifstream in(filename.c_str());
	if (!in.is_open())
	{
		std::cout << "Error in getNodeValues() - Could not find file for direct node values..." << std::endl;
		//return;
	}

	std::stringstream str_in;
	std::string line("");
	size_t idx(0);
	double val(0);

	while ( getline(in, line) )
	{
		if (line.find("#STOP") != std::string::npos)
			return;
		str_in << line;
		str_in >> idx >> val;
		node_values.push_back(std::pair<size_t, double>(idx, val));
		str_in.clear();
	}
}
