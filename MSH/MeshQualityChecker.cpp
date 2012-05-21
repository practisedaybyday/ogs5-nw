/*
 * MeshQualityChecker.cpp
 *
 *  Created on: Dec 8, 2010
 *      Author: TF
 */

#include "MeshQualityChecker.h"
#include "msh_elem.h"
#include <cmath>

namespace MeshLib
{
MeshQualityChecker::MeshQualityChecker(CFEMesh const* const mesh) :
	_min (), _max (), _mesh (mesh)
{
	if (_mesh)
		_mesh_quality_measure.resize ((_mesh->getElementVector()).size(), QualityType());
}

BASELIB::Histogram<double> MeshQualityChecker::getHistogram (size_t nclasses) const
{
	if (nclasses == 0) {
		// simple suggestion: number of classes with Sturges criterion
		nclasses = static_cast<size_t>(1 + 3.3 * log (static_cast<float>((_mesh->getElementVector()).size())));
	}

	return BASELIB::Histogram<double>(getMeshQuality(-1.0), nclasses, true);
}

void MeshQualityChecker::errorMsg (CElem* elem, size_t idx) const
{
	std::cout <<
	"Error in MeshQualityChecker::check() - Calculated value of element is below double precision minimum."
			  << std::endl;
	std::cout << "Points of " << MshElemType2String(elem->GetElementType()) << "-Element " <<
	idx << ": " << std::endl;
	for (int i(0); i < elem->GetVertexNumber(); i++)
		std::cout << "\t Node " << i << " " <<
		GEOLIB::Point((elem->GetNode(i))->getData()) << std::endl;
}

namespace detail
{
	struct QualityTypeToDoubleConverter
	{
		QualityTypeToDoubleConverter(const double& no_quality_value)
			: _no_quality_value(no_quality_value) { }
		double
		operator()(const boost::optional<double>& opt_value)
		{
			return (opt_value ? *opt_value : _no_quality_value);
		}

		private:
		const double _no_quality_value;
	};
}   // namespace detail

std::vector<double>
MeshQualityChecker::getMeshQuality (const double no_quality_value) const
{
	std::vector<double> q;
	q.reserve(_mesh_quality_measure.size());
	std::transform(
		_mesh_quality_measure.begin(),
		_mesh_quality_measure.end(),
		std::back_inserter(q),
		detail::QualityTypeToDoubleConverter(no_quality_value));
	return q;
}

double MeshQualityChecker::getMinValue(double no_quality_value) const
{
	return (_min ? *_min : no_quality_value);
}

double MeshQualityChecker::getMaxValue(double no_quality_value) const
{
	return (_max ? *_max : no_quality_value);
}

} // end namespace MeshLib
