/*
 * Polygon.h
 *
 *  Created on: Jun 21, 2010
 *      Author: TF
 */

#ifndef POLYGON_H_
#define POLYGON_H_

// STL
#include <list>

#include "Polyline.h"

namespace GEOLIB {

/**
 * \ingroup GEOLIB
 */

class Polygon : public Polyline
{
public:
	Polygon(const Polyline &ply);
	virtual ~Polygon();

	/**
	 * checks if the given point is inside the polygon
	 * @param pnt the Point
	 * @return if point is inside the polygon true, else false
	 */
	bool isPntInPolygon (const GEOLIB::Point& pnt) const;
	/**
	 * wrapper for method isPntInPolygon (const GEOLIB::Point&)
	 * @param x x coordinate of point
	 * @param y y coordinate of point
	 * @param z z coordinate of point
	 * @return if point is inside the polygon true, else false
	 */
	bool isPntInPolygon (double x, double y, double z) const;
	bool isPolylineInPolygon (const Polyline& ply) const;
	GEOLIB::Point* getIntersectionPointPolygonLine (GEOLIB::Point const & a, GEOLIB::Point const & b) const;
	void computeListOfSimplePolygons ();
	const std::list<Polygon*>& getListOfSimplePolygons ();

private:
	void splitPolygonAtIntersection (std::list<Polygon*>::iterator polygon_it);
	void splitPolygonAtPoint (std::list<Polygon*>::iterator polygon_it);
	double _maxx;
	double _maxy;
	std::list<Polygon*> _simple_polygon_list;
};

}

#endif /* POLYGON_H_ */
