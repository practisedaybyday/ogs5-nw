/**
 * \file BoundaryCondition.cpp
 * 2011/08/30 KR inital implementation
 *
 */

#include "BoundaryCondition.h"
#include "rf_bc_new.h"


BoundaryCondition::BoundaryCondition(const CBoundaryCondition &bc, const std::string &geometry_name)
: FEMCondition(geometry_name, bc.getProcessType(), bc.getProcessPrimaryVariable(), bc.getGeoType(), bc.getGeoName(),
			   bc.getProcessDistributionType(), FEMCondition::BOUNDARY_CONDITION)
{
	if (this->getProcessDistributionType() == FiniteElement::CONSTANT || this->getProcessDistributionType() == FiniteElement::CONSTANT_NEUMANN)
		this->setDisValue(bc.getGeoNodeValue());
	else if (this->getProcessDistributionType() == FiniteElement::LINEAR || this->getProcessDistributionType() == FiniteElement::LINEAR_NEUMANN)
		this->setLinearDisValues(bc.getPointsWithDistribedBC(), bc.getDistribedBC());
}
