
#include "DistributionTools.h"

#include "InterpolationAlgorithms/PiecewiseLinearInterpolation.h"
#include "mathlib.h"
#include "msh_mesh.h"
#include "LinearFunctionData.h"

void setDistributionConstant(MeshLib::CFEMesh &/*msh*/, std::vector<long> &vec_node_ids, std::vector<double> &vec_node_values, double val)
{
	const size_t nodes_vector_length = vec_node_ids.size();
	for (size_t i = 0; i < nodes_vector_length; i++)
		vec_node_values[i] = val;
}

void setDistributionGradient(MeshLib::CFEMesh &msh, std::vector<long> &vec_node_ids, std::vector<double> &vec_node_values, double ref_depth, double ref_value, double gradient)
{
	const size_t nodes_vector_length = vec_node_ids.size();
	for (size_t i = 0; i < nodes_vector_length; i++) {
		long node_id = vec_node_ids[i];
		vec_node_values[i] = ref_depth * (ref_depth - msh.nod_vector[node_id]->getData()[2]) + ref_value;
	}
}

void setDistributionFunction(MeshLib::CFEMesh &msh, std::vector<long> &vec_node_ids, std::vector<double> &vec_node_values, LinearFunctionData &linear_f)
{
	const size_t nodes_vector_length = vec_node_ids.size();
	for (size_t i = 0; i < nodes_vector_length; i++) {
		long node_id = vec_node_ids[i];
		double const* const coords(msh.nod_vector[node_id]->getData());
		vec_node_values[i] = linear_f.getValue(coords[0], coords[1], coords[2]);
	}
}

// Interpolation of polygon values to nodes_on_sfc
void SurfaceInterpolation(Surface* m_surface, MeshLib::CFEMesh* m_msh, std::vector<long>&nodes_on_sfc, std::vector<double>&node_value_vector)
{
	//----------------------------------------------------------------------
	const double Tol = m_msh->getMinEdgeLength();
	double gC[3], p1[3], p2[3], vn[3], unit[3], NTri[3];
	//
	std::vector<CGLPolyline*>::iterator p = m_surface->polyline_of_surface_vector.begin();

	for (size_t j=0; j<nodes_on_sfc.size(); j++)
	{
		double const* const pn (m_msh->nod_vector[nodes_on_sfc[j]]->getData());
		node_value_vector[j] = 0.0;
		bool Passed = false;
		// nodes close to first polyline
		p = m_surface->polyline_of_surface_vector.begin();
		while (p != m_surface->polyline_of_surface_vector.end())
		{
			CGLPolyline* m_polyline = *p;
			// Gravity center of this polygon
			for (int i = 0; i < 3; i++)
				gC[i] = 0.0;
			vn[2] = 0.0;
			size_t nPointsPly = m_polyline->point_vector.size();
			for (size_t i = 0; i < nPointsPly; i++)
			{
				gC[0] += m_polyline->point_vector[i]->x;
				gC[1] += m_polyline->point_vector[i]->y;
				gC[2] += m_polyline->point_vector[i]->z;
				vn[2] += m_polyline->point_vector[i]->getPropert();
			}
			for (size_t i = 0; i < 3; i++)
				gC[i] /= (double) nPointsPly;
			// BC value at center is an average of all point values of polygon
			vn[2] /= (double) nPointsPly;
			// Area of this polygon by the gravity center
			for (size_t i = 0; i < nPointsPly; i++)
			{
				p1[0] = m_polyline->point_vector[i]->x;
				p1[1] = m_polyline->point_vector[i]->y;
				p1[2] = m_polyline->point_vector[i]->z;
				size_t k = i + 1;
				if (i == nPointsPly - 1)
					k = 0;
				p2[0] = m_polyline->point_vector[k]->x;
				p2[1] = m_polyline->point_vector[k]->y;
				p2[2] = m_polyline->point_vector[k]->z;
				vn[0] = m_polyline->point_vector[i]->getPropert();
				vn[1] = m_polyline->point_vector[k]->getPropert();

				double Area1 = fabs(ComputeDetTri(p1, gC, p2));

				// Check if pn is in the triangle by points (p1, gC, p2)
				double Area2 = fabs(ComputeDetTri(p2, gC, pn));
				unit[0] = fabs(ComputeDetTri(gC, p1, pn));
				unit[1] = fabs(ComputeDetTri(p1, p2, pn));
				Area2 += unit[0] + unit[1];
				if (fabs(Area1 - Area2) < Tol)
				{
					// Interpolation within a triangle (p1,p2,gC)
					// Shape function
					for (size_t l = 0; l < 2; l++)
						unit[l] /= Area1;
					ShapeFunctionTri(NTri, unit);
					for (size_t l = 0; l < 3; l++)
						node_value_vector[j] += vn[l] * NTri[l];
					Passed = true;
					break;
				}
			}
			//
			p++;
			if (Passed)
				break;
		}                         // while
	}                                     //j
}

void setDistributionLinearPolyline(MeshLib::CFEMesh &msh, std::vector<long> &vec_node_ids, std::vector<double> &vec_node_values, GEOLIB::Polyline const* ply, std::vector<double> &DistribedBC, std::vector<int> &PointsHaveDistribedBC)
{
	std::vector<double> nodes_as_interpol_points;
	msh.getPointsForInterpolationAlongPolyline (ply, nodes_as_interpol_points);
	double msh_min_edge_length = msh.getMinEdgeLength();
	msh.setMinEdgeLength(msh_min_edge_length);

	std::vector<double> interpolation_points;
	std::vector<double> interpolation_values;
	for (size_t i(0); i < DistribedBC.size(); i++)
	{
		for (size_t j = 0; j < ply->getNumberOfPoints(); j++)
		{
			if (PointsHaveDistribedBC[i] == (int)ply->getPointID(j))
			{
				if (std::abs(DistribedBC[i]) < std::numeric_limits<double>::epsilon())
					DistribedBC[i] = 1.0e-20;
				interpolation_points.push_back (ply->getLength(j));
				interpolation_values.push_back (DistribedBC[i]);
				break;
			}
		}
	}
	MathLib::PiecewiseLinearInterpolation (
	        interpolation_points,
	        interpolation_values,
	        nodes_as_interpol_points,
	        vec_node_values);
}

void setDistributionLinearSurface(MeshLib::CFEMesh &msh, std::vector<long> &vec_node_ids, std::vector<double> &vec_node_values, Surface* m_surface, std::vector<double> &DistribedBC, std::vector<int> &PointsHaveDistribedBC)
{
	std::vector<CGLPolyline*>::iterator p = m_surface->polyline_of_surface_vector.begin();
	p = m_surface->polyline_of_surface_vector.begin();
	while (p != m_surface->polyline_of_surface_vector.end()) {
		CGLPolyline* m_polyline = *p;
		for (size_t i(0); i < DistribedBC.size(); i++) {
			for (size_t j = 0; j < m_polyline->point_vector.size(); j++)
				if (PointsHaveDistribedBC[i] == m_polyline->point_vector[j]->id) {
					if (std::abs(DistribedBC[i]) < std::numeric_limits<double>::epsilon())
						DistribedBC[i] = 1.0e-20;
					m_polyline-> point_vector[j]->setPropert(DistribedBC[i]);
					break;
				}
		}
		p++;
	}
	SurfaceInterpolation(m_surface, &msh, vec_node_ids, vec_node_values);

}

void setDistribution(DistributionData &dis_data, MeshLib::CFEMesh &msh, std::vector<long> &vec_node_ids, std::vector<double> &vec_node_values)
{
	if (dis_data.dis_type == FiniteElement::CONSTANT) {
		setDistributionConstant(msh, vec_node_ids, vec_node_values, dis_data.dis_parameters[0]);
	} else if (dis_data.dis_type == FiniteElement::LINEAR) {
		if (dis_data.geo_type==GEOLIB::POLYLINE) {
			GEOLIB::Polyline const* ply(static_cast<const GEOLIB::Polyline*> (dis_data.geo_obj));
			setDistributionLinearPolyline(msh, vec_node_ids, vec_node_values, ply, dis_data._DistribedBC, dis_data._PointsHaveDistribedBC);
		} else if (dis_data.geo_type==GEOLIB::SURFACE) {
			Surface* m_surface = dis_data.surface;
			setDistributionLinearSurface(msh, vec_node_ids, vec_node_values, m_surface, dis_data._DistribedBC, dis_data._PointsHaveDistribedBC);
		}
	} else if (dis_data.dis_type == FiniteElement::GRADIENT) {
		setDistributionGradient(msh, vec_node_ids, vec_node_values, dis_data.dis_parameters[0], dis_data.dis_parameters[1], dis_data.dis_parameters[2]);
	} else if (dis_data.dis_type == FiniteElement::FUNCTION) {
		setDistributionFunction(msh, vec_node_ids, vec_node_values, *dis_data.linear_f);
	}

}
