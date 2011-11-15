/**
 * \file FEMCondition.h
 * 25/11/2010 KR inital implementation
 *
 */

#ifndef FEMCONDITION_H
#define FEMCONDITION_H

#include "DistributionInfo.h"
#include "GeoInfo.h"
#include "GeoObject.h"
#include "Point.h"
#include "ProcessInfo.h"

#include <vector>

class CBoundaryCondition;
class CInitialCondition;
class CSourceTerm;

//class GEOObjects;
class GridAdapter;
#include "GEOObjects.h"

/**
 * \brief Adapter class for handling FEM Conditions in the user Interface
 */
class FEMCondition : public ProcessInfo, public GeoInfo, public DistributionInfo
{
public:
	/// Specifier for types of FEM Conditions
	enum CondType {
		UNSPECIFIED        = 0,
		BOUNDARY_CONDITION = 1,
		INITIAL_CONDITION  = 2,
		SOURCE_TERM        = 3
	};

	FEMCondition(const std::string &geometry_name, CondType t = UNSPECIFIED);
	FEMCondition(const std::string &geometry_name,
					FiniteElement::ProcessType pt = FiniteElement::INVALID_PROCESS,
					FiniteElement::PrimaryVariable pv =	FiniteElement::INVALID_PV,
					GEOLIB::GEOTYPE gt = GEOLIB::INVALID,
					const std::string &gn = "[unspecified]",
					FiniteElement::DistributionType dt = FiniteElement::INVALID_DIS_TYPE, CondType ct = UNSPECIFIED);
	FEMCondition(const FEMCondition &cond, CondType t);

	~FEMCondition() {}

	/// Returns the type of the FEM Condition (i.e. BC, IC or ST)
	CondType getCondType() const { return _type; }

	/// Returns the value(s) for the distribution
	const std::vector<double> getDisValue() const { return _disValue; }

	/// Returns the name of the geo-object the condition is assigned to. This object is part of the associated geometry.
	const std::string& getGeoName() const { return _geoName; }

	/// Returns the name of the associated geometry.
	const std::string& getAssociatedGeometryName() const { return _associated_geometry; }

	/// Sets a vector of values specifying the distribution.
	void setDisValue(const std::vector<double> &disValue)
	{
		for (size_t i = 0; i < disValue.size(); i++)
			_disValue.push_back(disValue[i]);
	}

	/// Sets a vector of values specifying the distribution.
	void setLinearDisValues(const std::vector<int> &point_ids,
	                        const std::vector<double> &point_values);

	/// Convenience method for setting a single value specifying the distribution.
	void setDisValue(double disValue) { _disValue.push_back(disValue); }

	/// Sets the name of the geo-object the condition is assigned to.
	void setGeoName(std::string geoName) { _geoName = geoName; }

	/// Returns the type of the FEM condition as a string.
	static std::string condTypeToString(CondType type);

protected:
	CondType _type;
	std::string _geoName;
	std::vector<double> _disValue;
	std::string _associated_geometry;
};

#endif //FEMCONDITION_H
