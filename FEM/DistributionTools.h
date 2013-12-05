
#ifndef DISTRIBUTIONTOOLS_H_
#define DISTRIBUTIONTOOLS_H_

#include <vector>
#include <string>

#include "GeoType.h"
#include "FEMEnums.h"

namespace GEOLIB
{
class GeoObject;
}
class Surface;

namespace MeshLib
{
class CFEMesh;
}

class LinearFunctionData;

struct DistributionData
{
	GEOLIB::GEOTYPE geo_type;
	const GEOLIB::GeoObject* geo_obj;
	Surface* surface;
	FiniteElement::DistributionType dis_type;
	std::vector<double> dis_parameters;
	LinearFunctionData* linear_f;
	std::vector<int> _PointsHaveDistribedBC;
	std::vector<double> _DistribedBC;
};

void setDistribution(DistributionData &dis_data, MeshLib::CFEMesh &msh, std::vector<long> &vec_node_ids, std::vector<double> &vec_node_values);

#endif                                            /* DISTRIBUTIONTOOLS_H_ */
