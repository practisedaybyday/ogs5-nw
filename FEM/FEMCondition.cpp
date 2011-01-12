/**
 * \file FEMCondition.cpp
 * 25/11/2010 KR inital implementation
 *
 */

#include "FEMCondition.h"


FEMCondition::FEMCondition(CondType t) 
: _type(t), _geoObject(NULL), _geoName("")
{
	this->setProcessType(INVALID_PROCESS);
	this->setProcessPrimaryVariable(INVALID_PV);
	this->setGeoType(GEOLIB::INVALID);
	this->setProcessDistributionType(FiniteElement::INVALID_DIS_TYPE);
}
