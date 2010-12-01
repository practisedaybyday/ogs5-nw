/*
 * \file AxisAlignedBoundingBox.h
 *
 *  Created on: April 22, 2010
 *      Author: TF
 */

#ifndef AXISALIGNEDBOUNDINGBOX_H_
#define AXISALIGNEDBOUNDINGBOX_H_

namespace GEOLIB {

/**
 *
 * \ingroup GEOLIB
 *
 * \brief Class AABB is a bounding box around a given geometric entity
 * */
class AABB
{
public:
	/**
	 * construction of object, initialization the axis aligned bounding box
	 * */
	AABB ();

	/**
	 * update axis aligned bounding box
	 */
	void update (double x, double y, double z);

	/**
	 * update axis aligned bounding box
	 */
	void update (const double *pnt)
	{
		update (pnt[0], pnt[1], pnt[2]);
	}

	/**
	 * wrapper for GEOLIB::Point
	 */
	bool containsPoint (const double *pnt) const;

	/**
	 * check if point described by its coordinates x, y, z is in
	 * the axis aligned bounding box
	 */
	bool containsPoint (double x, double y, double z) const;


private:
	double _min[3];
	double _max[3];
};

} // end namespace

#endif /* AXISALIGNEDBOUNDINGBOX_H_ */
