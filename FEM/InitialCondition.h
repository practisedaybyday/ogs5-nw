/**
 * \file InitialCondition.h
 * 2011/08/30 KR inital implementation
 *
 */

#ifndef INITIALCONDITION_H
#define INITIALCONDITION_H

#include "FEMCondition.h"


class InitialCondition : public FEMCondition
{
public:
	InitialCondition(const std::string &geometry_name) : FEMCondition(geometry_name, FEMCondition::INITIAL_CONDITION) {};
	InitialCondition(const CInitialCondition &ic, const std::string &geometry_name);
	~InitialCondition() {};
};

#endif //INITIALCONDITION_H
