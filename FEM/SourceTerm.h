/**
 * \file SourceTerm.h
 * 2011/08/30 KR inital implementation
 *
 */

#ifndef SOURCETERM_H
#define SOURCETERM_H

#include "FEMCondition.h"

class SourceTerm : public FEMCondition
{
public:
	SourceTerm(const std::string &geometry_name)
		: FEMCondition(geometry_name, FEMCondition::SOURCE_TERM), _tim_type(0) {}
	SourceTerm(const CSourceTerm &st, const std::string &geometry_name);
	~SourceTerm() {}

	size_t getTimType() const {return _tim_type; }
	void setTimType(size_t value) { _tim_type = value; }

	static std::vector<FEMCondition*> createDirectSourceTerms(
	        const std::vector<CSourceTerm*> &st_vector,
	        const std::string &geo_name);

private:
	static void getDirectNodeValues(const std::string &filename,
	                                std::vector< std::pair<size_t, double> > &nodes_values);
	size_t _tim_type;
};

#endif //SOURCETERM_H
