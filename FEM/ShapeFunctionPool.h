/*! \file ShapeFunctionPool.h
     \brief Compute shape functions and their gradients with respect to
	  the local coodinates, and store the results

	 \author Wenqing Wang
	 \date Feb. 2015
*/

#ifndef OGS_SHAPEFUNCTIONPOOL_H
#define OGS_SHAPEFUNCTIONPOOL_H

#include <vector>

#include "MSHEnums.h"

namespace FiniteElement
{
class CElement;

class ShapeFunctionPool
{
public:
	ShapeFunctionPool(const std::vector<MshElemType::type>& elem_types,
		              CElement& quadrature);
	~ShapeFunctionPool();

	/// Get shape function values of an element type
	double* getShapeFunctionValues(const MshElemType::type elem_type) const;

	/// Get the values of the gradient of shape function of an element type
	double* getGradShapeFunctionValues(const MshElemType::type elem_type) const;

private:
	/// Results of shape functions of all integration points.
	std::vector<double*> _shape_fucntion; 
	/// Results of the gradient of shape functions with respect to
	/// local coordinates of all integration points.
	std::vector<double*> _grad_shape_fucntion; 

	void computeQuadratures(const std::vector<MshElemType::type>& elem_types,
		                    const int num_elem_nodes[2][MshElemType::LAST],
							const int dim_elem[],
		                    CElement& quadrature);
};
} // end namespace

#endif
