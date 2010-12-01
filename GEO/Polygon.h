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

	bool isPntInPolygon (const GEOLIB::Point& pnt) const;
	bool isPolylineInPolygon (const Polyline& ply) const;
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
