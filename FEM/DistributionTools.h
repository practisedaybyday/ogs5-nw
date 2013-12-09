
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
	std::string geo_name;
	const GEOLIB::GeoObject* geo_obj;
	FiniteElement::DistributionType dis_type;
	std::vector<double> dis_parameters;
	LinearFunctionData* linear_f;
	std::vector<int> _PointsHaveDistribedBC;
	std::vector<double> _DistribedBC;
	std::string mesh_type_name;
	size_t mesh_node_id;
	int mat_id;

	DistributionData()
	: geo_type(GEOLIB::INVALID), geo_obj(NULL), dis_type(FiniteElement::INVALID_DIS_TYPE), linear_f(NULL),
	  mesh_node_id(-1), mat_id(-1)
	{}
};

void getNodesOnDistribution(DistributionData &dis_data, MeshLib::CFEMesh &msh, std::vector<long> &vec_node_ids);

void setDistribution(DistributionData &dis_data, MeshLib::CFEMesh &msh, std::vector<long> &vec_node_ids, std::vector<double> &vec_node_values);

#endif                                            /* DISTRIBUTIONTOOLS_H_ */
