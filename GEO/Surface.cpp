/*
 * Surface.cpp
 *
 *  Created on: Apr 22, 2010
 *      Author: TF
 */

#include "Surface.h"
#include "AxisAlignedBoundingBox.h"

namespace GEOLIB {

Surface::Surface (const std::vector<Point*> &pnt_vec) :
	GeoObject(), m_sfc_pnts(pnt_vec), bv (new AABB())
{}

Surface::~Surface ()
{
	for (size_t k(0); k<m_sfc_triangles.size(); k++)
		delete m_sfc_triangles[k];
	delete bv;
}

void Surface::addTriangle (size_t pnt_a, size_t pnt_b, size_t pnt_c)
{
	assert (pnt_a < m_sfc_pnts.size() && pnt_b < m_sfc_pnts.size() && pnt_c < m_sfc_pnts.size());
	m_sfc_triangles.push_back (new Triangle(m_sfc_pnts, pnt_a, pnt_b, pnt_c));
	bv->update (*m_sfc_pnts[pnt_a]);
	bv->update (*m_sfc_pnts[pnt_b]);
	bv->update (*m_sfc_pnts[pnt_c]);
}

size_t Surface::getNTriangles () const
{
	return m_sfc_triangles.size();
}

const Triangle* Surface::operator[] (size_t i) const
{
	assert (i < m_sfc_triangles.size());
	return m_sfc_triangles[i];
}

bool Surface::isPntInBV (const double *pnt) const
{
	return bv->containsPoint (pnt);
}

bool Surface::isPntInSfc (const double *pnt) const
{
	bool nfound (true);
	for (size_t k(0); k<m_sfc_triangles.size() && nfound; k++) {
		if (m_sfc_triangles[k]->containsPoint (pnt)) nfound = false;
	}
	return !nfound;
}

} // end namespace
