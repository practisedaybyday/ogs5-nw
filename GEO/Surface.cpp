/*
 * Surface.cpp
 *
 *  Created on: Apr 22, 2010
 *      Author: TF
 */

#include <list>

#include "display.h"

// GEOLIB
#include "AxisAlignedBoundingBox.h"
#include "Polygon.h"
#include "Surface.h"

// MathLib
#include "AnalyticalGeometry.h"
#include "EarClippingTriangulation.h"

namespace GEOLIB
{
Surface::Surface (const std::vector<Point*> &pnt_vec) :
	GeoObject(), _sfc_pnts(pnt_vec), _bv(), _sfc_grid(NULL)
{}

Surface::~Surface ()
{
	delete _sfc_grid;
	for (size_t k(0); k < _sfc_triangles.size(); k++)
		delete _sfc_triangles[k];
}

void Surface::addTriangle (size_t pnt_a, size_t pnt_b, size_t pnt_c)
{
	assert (pnt_a < _sfc_pnts.size() && pnt_b < _sfc_pnts.size() && pnt_c < _sfc_pnts.size());
	_sfc_triangles.push_back (new Triangle(_sfc_pnts, pnt_a, pnt_b, pnt_c));
	_bv.update (*_sfc_pnts[pnt_a]);
	_bv.update (*_sfc_pnts[pnt_b]);
	_bv.update (*_sfc_pnts[pnt_c]);
	delete _sfc_grid;
	_sfc_grid=NULL;
}

Surface* Surface::createSurface(const Polyline &ply)
{
	if (!ply.isClosed())
	{
		std::cout << "Error in Surface::createSurface() - Polyline is not closed..." << std::endl;
		return NULL;
	}

	if (ply.getNumberOfPoints() > 2)
	{
		// create empty surface
		Surface* sfc(new Surface(ply.getPointsVec()));

		Polygon* polygon (new Polygon (ply));
		polygon->computeListOfSimplePolygons ();

		// create surfaces from simple polygons
		const std::list<GEOLIB::Polygon*>& list_of_simple_polygons (
		        polygon->getListOfSimplePolygons());
		for (std::list<GEOLIB::Polygon*>::const_iterator simple_polygon_it (
		             list_of_simple_polygons.begin());
		     simple_polygon_it != list_of_simple_polygons.end(); ++simple_polygon_it)
		{
			std::list<GEOLIB::Triangle> triangles;
			ScreenMessage("-> triangulation of surface: ... \n");
			MathLib::EarClippingTriangulation(*simple_polygon_it, triangles);
			ScreenMessage("-> done - %d triangles\n", triangles.size());

			// add Triangles to Surface
			std::list<GEOLIB::Triangle>::const_iterator it (triangles.begin());
			while (it != triangles.end())
			{
				sfc->addTriangle ((*it)[0], (*it)[1], (*it)[2]);
				it++;
			}
		}
		delete polygon;
		return sfc;
	}
	else
	{
		std::cout <<
		"Error in Surface::createSurface() - Polyline consists of less than three points and therefore cannot be triangulated..."
		          << std::endl;
		return NULL;
	}
}

size_t Surface::getNTriangles () const
{
	return _sfc_triangles.size();
}

const Triangle* Surface::operator[] (size_t i) const
{
	assert (i < _sfc_triangles.size());
	return _sfc_triangles[i];
}

bool Surface::isPntInBV (const double* pnt, double eps) const
{
	return _bv.containsPoint (pnt, eps);
}

void Surface::initSurfaceGrid()
{
	if (_sfc_grid == NULL) {
		_sfc_grid = new SurfaceGrid(this);
	}
}

bool Surface::isPntInSfc (const double* pnt, double eps) const
{
	return _sfc_grid->isPntInSurface(pnt, eps);
}

void Surface::removeTriangles(const GEOLIB::AABB &bb)
{
	// find triangles located inside of the given bounding box
	std::vector<Triangle*> new_tris;
	for (size_t i(0); i < _sfc_triangles.size(); i++)
	{
		Triangle* tri = _sfc_triangles[i];
		bool intersecting = false;
		for (size_t k(0); k < 3; k++)
		{
			if (bb.containsPoint(*tri->getPoint(k))) {
				intersecting = true;
				break;
			}
		}
		if (intersecting) {
			new_tris.push_back(tri);
		}
		else
		{
			for (unsigned j=0; j<3; j++) {
				delete _sfc_pnts[(*tri)[j]];
				(*const_cast<std::vector<Point*>*>(&_sfc_pnts))[(*tri)[j]] = NULL;
			}
			delete tri;
		}
	}

	if (new_tris.size()==_sfc_triangles.size()) {
		return;
	}

	// initialize surface data
	_sfc_triangles.clear();
	_bv.clear();
	delete _sfc_grid;
	_sfc_grid=NULL;

	// update
	for (size_t i(0); i < new_tris.size(); i++)
	{
		_sfc_triangles.push_back(new_tris[i]);
		_bv.update (*_sfc_pnts[(*new_tris[i])[0]]);
		_bv.update (*_sfc_pnts[(*new_tris[i])[1]]);
		_bv.update (*_sfc_pnts[(*new_tris[i])[2]]);
	}

}

#if 0
void Surface::removeTriangles()
{
	for (size_t k(0); k < _sfc_triangles.size(); k++)
		delete _sfc_triangles[k];
	_sfc_triangles.clear();
}
#endif

} // end namespace
