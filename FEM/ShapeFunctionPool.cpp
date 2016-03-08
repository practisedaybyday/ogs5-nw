/*! \file ShapeFunctionPool.cpp
     \brief Compute shape functions and their gradients with respect to
	  the local coodinates, and store the results

	 \author Wenqing Wang
	 \date Feb. 2015
*/

#include "ShapeFunctionPool.h"

#include <cassert>     /* assert */
#include "fem_ele.h"

namespace FiniteElement
{
ShapeFunctionPool::ShapeFunctionPool(
			const std::vector<MshElemType::type>& elem_types,
			CElement& quadrature, const int num_sample_gs_pnts)
{
	int num_elem_nodes[2][MshElemType::LAST];
	int dim_elem[MshElemType::LAST];

	const std::size_t n_ele_types = elem_types.size();
	_shape_fucntion.reserve(n_ele_types);
	_grad_shape_fucntion.reserve(n_ele_types);

	for (std::size_t i=0; i<n_ele_types; i++)
	{
		_shape_fucntion.push_back(NULL);
	}

	for (std::size_t i=0; i<n_ele_types; i++)
	{
		_grad_shape_fucntion.push_back(NULL);
	}

	int id = MshElemType::LINE - 1;
	num_elem_nodes[0][id] = 2;
	num_elem_nodes[1][id] = 3;
	dim_elem[id] = 1;

	id = MshElemType::QUAD - 1;
	num_elem_nodes[0][id] = 4;
	num_elem_nodes[1][id] = 9;
	dim_elem[id] = 2;

	id = MshElemType::QUAD8 - 1;
	num_elem_nodes[0][id] = 4;
	num_elem_nodes[1][id] = 8;
	dim_elem[id] = 2;

	id = MshElemType::TRIANGLE - 1;
	num_elem_nodes[0][id] = 3;
	num_elem_nodes[1][id] = 6;
	dim_elem[id] = 2;

	id = MshElemType::HEXAHEDRON - 1;
	num_elem_nodes[0][id] = 8;
	num_elem_nodes[1][id] = 20;
	dim_elem[id] = 3;

	id = MshElemType::TETRAHEDRON - 1;
	num_elem_nodes[0][id] = 4;
	num_elem_nodes[1][id] = 10;
	dim_elem[id] = 3;

	id = MshElemType::PRISM - 1;
	num_elem_nodes[0][id] = 6;
	num_elem_nodes[1][id] = 15;
	dim_elem[id] = 3;

	id = MshElemType::PYRAMID - 1;
	num_elem_nodes[0][id] = 5;
	num_elem_nodes[1][id] = 13;
	dim_elem[id] = 3;

	//std::vector<int> elem_type_ids
	for (std::size_t i=0; i<elem_types.size(); i++)
	{
		const int type_id = elem_types[i] - 1;
		if (elem_types[type_id] == MshElemType::INVALID)
			continue;
		// Set number of integration points.
		quadrature.SetGaussPointNumber(num_sample_gs_pnts);
		quadrature.ConfigShapefunction(elem_types[i]);

		int num_int_pnts = quadrature.GetNumGaussPoints();
		const int size_shape_fct = num_elem_nodes[quadrature.getOrder() - 1][type_id]
		                            * num_int_pnts;
		_shape_fucntion[type_id] = new double[size_shape_fct];
		_grad_shape_fucntion[type_id] = new double[dim_elem[type_id] * size_shape_fct];
	}

	computeQuadratures(elem_types, num_elem_nodes, dim_elem, quadrature, num_sample_gs_pnts);
}

ShapeFunctionPool::~ShapeFunctionPool()
{
	for (std::size_t i=0; i<_shape_fucntion.size(); i++)
	{
		if (_shape_fucntion[i])
			delete [] _shape_fucntion[i];
		_shape_fucntion[i] = NULL;
	}

	for (std::size_t i=0; i<_grad_shape_fucntion.size(); i++)
	{
		if (_grad_shape_fucntion[i])
			delete [] _grad_shape_fucntion[i];
		_grad_shape_fucntion[i] = NULL;
	}
}

void ShapeFunctionPool::
	computeQuadratures(const std::vector<MshElemType::type>& elem_types,
	                   const int num_elem_nodes[2][MshElemType::LAST],
					   const int dim_elem[],
		               CElement& quadrature, const int num_sample_gs_pnts)
{
	const int order = quadrature.getOrder();
	for (std::size_t i=0; i<elem_types.size(); i++)
	{
		const int type_id = elem_types[i] - 1;
		quadrature.ConfigShapefunction(elem_types[i]);

		if (elem_types[type_id] == MshElemType::INVALID)
			continue;

    	const int nnodes = num_elem_nodes[order-1][type_id];
    	const int elem_dim = dim_elem[type_id];
		double* shape_function_values = _shape_fucntion[type_id];
		double* dshape_function_values = _grad_shape_fucntion[type_id];

		// Set number of integration points.
		quadrature.SetGaussPointNumber(num_sample_gs_pnts);
		quadrature.ConfigShapefunction(elem_types[i]);
		for (int gp = 0; gp < quadrature.GetNumGaussPoints(); gp++)
	    {
			int gp_r, gp_s, gp_t;
			quadrature.SetGaussPoint(elem_types[i], gp, gp_r, gp_s, gp_t);
			double* shape_function_values_gs
				    = &shape_function_values[gp * nnodes];
			quadrature.ComputeShapefct(order, shape_function_values_gs);

			double* dshape_function_values_gs
				    = &dshape_function_values[gp * nnodes * elem_dim];
			quadrature.ComputeGradShapefctLocal(order, dshape_function_values_gs);
		}
	}
}

double* ShapeFunctionPool::
getShapeFunctionValues(const MshElemType::type elem_type) const
{
	assert(_shape_fucntion[elem_type-1]);
	return _shape_fucntion[elem_type-1];
}

double* ShapeFunctionPool::
getGradShapeFunctionValues(const MshElemType::type elem_type) const
{
	assert(_shape_fucntion[elem_type-1]);
	return _grad_shape_fucntion[elem_type-1];
}

} // end namespace


