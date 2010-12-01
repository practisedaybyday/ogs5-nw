/*
 * \file AnalyticalGeometry.cpp
 *
 *  Created on: Mar 17, 2010
 *      Author: TF
 */

#include <cmath>
#include <cstdlib> // for exit
#include <list>
#include <limits>
#include <fstream>

// Base
#include "swap.h"
#include "quicksort.h"

// GEO
#include "Polyline.h"
#include "Triangle.h"

// MathLib
#include "MathTools.h"
#include "AnalyticalGeometry.h"
#include "GaussAlgorithm.h"
#include "Matrix.h" // for transformation matrix
#include "max.h"

namespace MATHLIB {

size_t getOrientation (const double& p0_x, const double& p0_y,
	const double& p1_x, const double& p1_y,
	const double& p2_x, const double& p2_y)
{
	double h1 ((p1_x-p0_x)*(p2_y-p0_y));
	double h2 ((p2_x-p0_x)*(p1_y-p0_y));

	double tol (sqrt(std::numeric_limits<double>::min()));
	if (fabs (h1-h2) <= tol * max (fabs(h1), fabs(h2)))
		return COLLINEAR;
	if (h1-h2 > 0.0) return CCW;

	return CW;
}

size_t getOrientation (const GEOLIB::Point* p0, const GEOLIB::Point* p1, const GEOLIB::Point* p2)
{
	return getOrientation ((*p0)[0], (*p0)[1], (*p1)[0], (*p1)[1], (*p2)[0], (*p2)[1]);
}

bool lineSegmentIntersect (const GEOLIB::Point& a, const GEOLIB::Point& b,
		const GEOLIB::Point& c, const GEOLIB::Point& d,
		GEOLIB::Point& s)
{
	Matrix<double> mat(2,2);
	mat(0,0) = b[0] - a[0];
	mat(1,0) = b[1] - a[1];
	mat(0,1) = c[0] - d[0];
	mat(1,1) = c[1] - d[1];

	// check if vectors are parallel
	double eps (sqrt(std::numeric_limits<double>::min()));
	if (fabs(mat(1,1)) < eps) {
		// vector (D-C) is parallel to x-axis
		if (fabs(mat(0,1)) < eps) {
			// vector (B-A) is parallel to x-axis
			return false;
		}
	} else {
		// vector (D-C) is not parallel to x-axis
		if (fabs(mat(0,1)) >= eps) {
			// vector (B-A) is not parallel to x-axis
			// \$f(B-A)\f$ and \$f(D-C)\f$ are parallel iff there exists
			// a constant \f$c\f$ such that \$f(B-A) = c (D-C)\f$
			if (fabs (mat(0,0) / mat(0,1) - mat(1,0) / mat(1,1)) < eps * fabs (mat(0,0) / mat(0,1)))
				return false;
		}
	}

	double *rhs (new double[2]);
	rhs[0] = c[0] - a[0];
	rhs[1] = c[1] - a[1];

	GaussAlgorithm lu_solver (mat);
	lu_solver.execute (rhs);
	if (0 <= rhs[0] && rhs[0] <= 1.0 && 0 <= rhs[1] && rhs[1] <= 1.0) {
		s[0] = a[0] + rhs[0] * (b[0] - a[0]);
		s[1] = a[1] + rhs[0] * (b[1] - a[1]);
		s[2] = a[2] + rhs[0] * (b[2] - a[2]);
		delete [] rhs;
		return true;
	} else delete [] rhs;
	return false;
}

bool lineSegmentsIntersect (const GEOLIB::Polyline* ply, size_t &idx0, size_t &idx1, GEOLIB::Point& intersection_pnt)
{
	size_t n_segs (ply->getNumberOfPoints() - 1);
	/**
	 * computing the intersections of all possible pairs of line segments of the given polyline
	 * as follows:
	 * let the segment \f$s_1 = (A,B)\f$ defined by \f$k\f$-th and \f$k+1\f$-st point
	 * of the polyline and segment \f$s_2 = (C,B)\f$ defined by \f$j\f$-th and
	 * \f$j+1\f$-st point of the polyline, \f$j>k+1\f$
	 */
	for (size_t k(0); k<n_segs-2; k++) {
		for (size_t j(k+2); j<n_segs; j++) {
			if (k!=0 || j<n_segs-1) {
				if (lineSegmentIntersect (*(*ply)[k], *(*ply)[k+1], *(*ply)[j], *(*ply)[j+1], intersection_pnt)) {
					idx0 = k;
					idx1 = j;
					return true;
				}
			}
		}
	}
	return false;
}

bool isPointInTriangle (const double pp[3], const double pa[3], const double pb[3], const double pc[3])
{
	MATHLIB::Vector a(pa[0],pa[1],pa[2]), b(pb[0],pb[1],pb[2]), c(pc[0],pc[1],pc[2]);
	MATHLIB::Vector p (pp[0],pp[1],pp[2]);
	// move triangle
	a -= p;
	b -= p;
	c -= p;

	// compute normal vectors for triangles pab and pbc
	MATHLIB::Vector u (a.Cross (b));
	MATHLIB::Vector v (b.Cross (c));
	// make sure they are both pointing in the same direction
	if (u.Dot(v) < 0.0) return false;
	// compute normal vector for triangle pca
	MATHLIB::Vector w (c.Cross (a));
	// make sure it points in the same direction as the first two
	if (u.Dot(w) < 0.0) return false;

	return true;
}

bool isPointInTriangle (const GEOLIB::Point* p,
		const GEOLIB::Point* a, const GEOLIB::Point* b, const GEOLIB::Point* c)
{
	return isPointInTriangle (p->getData(), a->getData(), b->getData(), c->getData());
}


// NewellPlane from book Real-Time Collision detection p. 494
void getNewellPlane(const std::vector<GEOLIB::Point*>& pnts, Vector &plane_normal,
		double& d)
{
	d = 0;
	Vector centroid;
	size_t n_pnts (pnts.size());
	for (size_t i(n_pnts - 1), j(0); j < n_pnts; i = j, j++) {
		plane_normal[0] += ((*(pnts[i]))[1] - (*(pnts[j]))[1])
				* ((*(pnts[i]))[2] + (*(pnts[j]))[2]); // projection on yz
		plane_normal[1] += ((*(pnts[i]))[2] - (*(pnts[j]))[2])
				* ((*(pnts[i]))[0] + (*(pnts[j]))[0]); // projection on xz
		plane_normal[2] += ((*(pnts[i]))[0] - (*(pnts[j]))[0])
				* ((*(pnts[i]))[1] + (*(pnts[j]))[1]); // projection on xy

		centroid += *(pnts[j]);
	}

	plane_normal *= 1.0 / plane_normal.Length();
	d = centroid.Dot(plane_normal) / n_pnts;
}


void rotatePointsToXY(Vector &plane_normal,
		std::vector<GEOLIB::Point*> &pnts)
{
	// *** some frequently used terms ***
	// sqrt (v_1^2 + v_2^2)
	double h0(sqrt(plane_normal[0] * plane_normal[0] + plane_normal[1]
			* plane_normal[1]));
	// 1 / sqrt (v_1^2 + v_2^2)
	double h1(1 / h0);
	// 1 / sqrt (h0 + v_3^2)
	double h2(1.0 / sqrt(h0 + plane_normal[2] * plane_normal[2]));

	Matrix<double> rot_mat(3, 3);
	// calc rotation matrix
	rot_mat(0, 0) = plane_normal[2] * plane_normal[0] * h2 * h1;
	rot_mat(0, 1) = plane_normal[2] * plane_normal[1] * h2 * h1;
	rot_mat(0, 2) = - h0 * h2;
	rot_mat(1, 0) = -plane_normal[1] * h1;
	rot_mat(1, 1) = plane_normal[0] * h1;;
	rot_mat(1, 2) = 0.0;
	rot_mat(2, 0) = plane_normal[0] * h2;
	rot_mat(2, 1) = plane_normal[1] * h2;
	rot_mat(2, 2) = plane_normal[2] * h2;

	double *tmp (new double[3]);
	size_t n_pnts(pnts.size());
	for (size_t k(0); k < n_pnts; k++) {
		tmp = rot_mat * pnts[k]->getData();
		for (size_t j(0); j < 3; j++)
			(*(pnts[k]))[j] = tmp[j];
	}

	tmp = rot_mat * plane_normal.getData();
	for (size_t j(0); j < 3; j++) plane_normal[j] = tmp[j];

	delete [] tmp;
}

bool isEar(const std::vector<GEOLIB::Point*>& pnts,
		const std::list<size_t>& vertex_list, size_t v0, size_t v1, size_t v2)
{
	for (std::list<size_t>::const_iterator it (vertex_list.begin ());
		it != vertex_list.end(); ++it) {
		if (*it != v0 && *it != v1 && *it != v2) {
			if (isPointInTriangle (pnts[*it], pnts[v0], pnts[v1], pnts[v2])) {
				return false;
			}
		}
	}
	return true;
}

void uniqueListInsert (std::list<size_t>& list, size_t element)
{
	// search element
	std::list<size_t>::const_iterator it;
	for (it = list.begin (); it != list.end(); it++) {
		if (*it == element) return;
	}
	// element not found -> insert
	list.push_back (element);
}

void earClippingTriangulation (std::vector<GEOLIB::Point*>& pnts, std::list<GEOLIB::Triangle> &triangles)
{
	size_t n_pnts (pnts.size()), orientation;
	std::list<size_t> vertex_list;
	for (size_t k(0); k<n_pnts; k++) vertex_list.push_back (k);

	// initialize lists
	std::list<size_t> convex_vertex_list, ear_list;
	// go through points checking ccw, cw or collinear order and identifying ears
	std::list<size_t>::iterator it (vertex_list.begin()), prev(vertex_list.end()), next;
	prev--;
	next = it;
	next++;
	while (next != vertex_list.end()) {
		orientation  = getOrientation (pnts[*prev], pnts[*it], pnts[*next]);
		if (orientation == COLLINEAR) {
			it = vertex_list.erase (it);
			it = next;
			next++;
		} else {
			if (orientation == CCW) {
				convex_vertex_list.push_back (*it);
				bool is_ear (true);
				for (std::list<size_t>::const_iterator my_it (vertex_list.begin());
					my_it != vertex_list.end() && is_ear; ++my_it) {
					if (my_it != it && my_it != prev && my_it != next) {
						if (isPointInTriangle (pnts[*my_it], pnts[*prev], pnts[*it], pnts[*next]))
							is_ear = false;
					}
				}
				if (is_ear) {
					ear_list.push_back (*it);
				}
			}
			prev = it;
			it = next;
			next++;
		}
	}

	next = vertex_list.begin();
	orientation = getOrientation (pnts[*prev], pnts[*it], pnts[*next]);
	if (orientation == COLLINEAR) {
		it = vertex_list.erase (it);
	}
	if (orientation == CCW) {
		convex_vertex_list.push_back (*it);
		bool is_ear (true);
		for (std::list<size_t>::const_iterator my_it (vertex_list.begin());
			my_it != vertex_list.end() && is_ear; ++my_it) {
			if (my_it != it && my_it != prev && my_it != next) {
				if (isPointInTriangle (pnts[*my_it], pnts[*prev], pnts[*it], pnts[*next]))
					is_ear = false;
			}
		}
		if (is_ear) {
			ear_list.push_back (*it);
		}
	}

	// *** clip an ear
	while (vertex_list.size() > 3) {
		// pop ear from list
		size_t ear = ear_list.front ();
		ear_list.pop_front();

		// remove ear from vertex_list
		bool nfound (true);
		it = vertex_list.begin();
		prev = vertex_list.end();
		prev--;
		while (nfound && it != vertex_list.end()) {
			if (*it == ear) {
				nfound = false;
				it = vertex_list.erase (it); // remove ear tip
				next = it;
				if (next == vertex_list.end()) {
					next = vertex_list.begin();
					prev = vertex_list.end();
					prev--;
				}
				// check the orientation of prevprev, prev, next
				std::list<size_t>::iterator prevprev;
				if (prev == vertex_list.begin ()) {
					prevprev = vertex_list.end();
				} else {
					prevprev = prev;
				}
				prevprev--;

				// apply changes to convex_vertex_list and ear_list looking "backward"
				orientation = getOrientation (pnts[*prevprev], pnts[*prev], pnts[*next]);
				if (orientation == CCW) {
					// prev is convex
					if (isEar (pnts, vertex_list, *prevprev, *prev, *next)) {
						// prev is an ear tip
						uniqueListInsert (convex_vertex_list, *prev);
						uniqueListInsert (ear_list, *prev);
					} else {
						uniqueListInsert (convex_vertex_list, *prev);
						// if necessary remove prev
						ear_list.remove (*prev);
					}
				} else {
					// prev is not convex => reflex or collinear
					convex_vertex_list.remove (*prev);
					ear_list.remove (*prev);
					if (orientation == COLLINEAR) {
						prev = vertex_list.erase (prev);
						prev--;
					}
				}

				// check the orientation of prev, next, nextnext
				std::list<size_t>::iterator nextnext, help_it (vertex_list.end ());
				help_it--;
				if (next == help_it) {
					nextnext = vertex_list.begin();
				} else {
					nextnext = next;
					nextnext++;
				}

				// apply changes to convex_vertex_list and ear_list looking "forward"
				orientation = getOrientation (pnts[*prev], pnts[*next], pnts[*nextnext]);
				if (orientation == CCW) {
					// next is convex
					if (isEar (pnts, vertex_list, *prev, *next, *nextnext)) {
						// next is an ear tip
						uniqueListInsert (convex_vertex_list, *next);
						uniqueListInsert (ear_list, *next);
					} else {
						uniqueListInsert (convex_vertex_list, *next);
						// if necessary remove *prev
						ear_list.remove (*next);
					}
				} else {
					// next is not convex => reflex or collinear
					convex_vertex_list.remove (*next);
					ear_list.remove (*next);
					if (orientation == COLLINEAR) {
						next = vertex_list.erase (next);
						if (next == vertex_list.end())
							next = vertex_list.begin();
					}
				}

			} else {
				prev = it;
				it++;
			}
		}
		// remove ear tip from convex_vertex_list
		convex_vertex_list.remove (ear);
		// add triangle
		triangles.push_back (GEOLIB::Triangle (pnts, *prev, ear, *next));
	}
	// add last triangle
	next = vertex_list.begin();
	prev = next;
	next++;
	it = next;
	next++;
	if (getOrientation (pnts[*prev], pnts[*it], pnts[*next]) == CCW)
		triangles.push_back (GEOLIB::Triangle (pnts, *prev, *it, *next));
	else
		triangles.push_back (GEOLIB::Triangle (pnts, *prev, *next, *it));
}

void earClippingTriangulationOfPolygon(const GEOLIB::Polyline* ply, std::list<GEOLIB::Triangle> &triangles)
{
	if (!ply->isClosed()) {
		std::cerr << "Polyline is not a Polygon!" << std::endl;
		return;
	}
	// copy points - last point is identical to the first
	size_t n_pnts (ply->getNumberOfPoints()-1);

	std::vector<GEOLIB::Point*> pnts;
	pnts.reserve (n_pnts);
	for (size_t k(0); k < n_pnts; k++)
		pnts.push_back (new GEOLIB::Point (*(ply->getPoint(k))));

	// calculate supporting plane
	Vector plane_normal;
	double d;
	// compute the plane normal
	getNewellPlane(pnts, plane_normal, d);

	double tol (sqrt(std::numeric_limits<double>::min()));
	if (fabs(plane_normal[0]) > tol || fabs(plane_normal[1]) > tol) {
		// rotate copied points into x-y-plane
		rotatePointsToXY(plane_normal, pnts);
	}

	// recompute the plane normal
	getNewellPlane(pnts, plane_normal, d);

	// check orientation and if CW change orientation
	size_t orientation (CCW);
	double z_axis[3] = {0.0, 0.0, 1.0};
	if (MATHLIB::scpr (plane_normal.getData(), z_axis, 3) <= 0.0) {
		// CW orientation of polygon
		orientation = CW;
		for (size_t k(0); k<n_pnts/2; k++) {
			BASELIB::swap (pnts[k], pnts[n_pnts-1-k]);
		}
	}

	std::list<GEOLIB::Triangle> temp_triangles;
	earClippingTriangulation (pnts, temp_triangles);

	const std::vector<GEOLIB::Point*>& ref_pnts_vec (ply->getPointsVec());
	std::list<GEOLIB::Triangle>::const_iterator it (temp_triangles.begin());
	while (it != temp_triangles.end()) {
		const size_t i0 (ply->getPointID ((*it)[0]));
		const size_t i1 (ply->getPointID ((*it)[1]));
		const size_t i2 (ply->getPointID ((*it)[2]));

		if (orientation == CCW) {
			triangles.push_back (GEOLIB::Triangle (ref_pnts_vec, i0, i1, i2));
		} else {
			triangles.push_back (GEOLIB::Triangle (ref_pnts_vec,
					ply->getPointID (n_pnts-1-((*it)[0])),
					ply->getPointID (n_pnts-1-((*it)[1])),
					ply->getPointID (n_pnts-1-((*it)[2]))));
		}
		it++;
	}

	// delete points
	for (size_t k(0); k < n_pnts; k++) delete pnts[k];
}

//#include "LinkedTriangle.h" // used by Delaunay-Triangulation
//void getCircumscribedSphereOfTriangle(const double a[3], const double b[3],
//		const double c[3], double middle_pnt[3], double& sqr_radius)
//{
//	// describe first plane E_0 in normal form (between a, b)
//	double n0[3], d0(0.0);
//	for (size_t k(0); k < 3; k++) {
//		n0[k] = b[k] - a[k];
//		d0 -= n0[k] * (a[k] + 0.5 * n0[k]);
//	}
//
//	// describe second plane E_1 in normal form (between a, c)
//	double n1[3], d1(0.0);
//	for (size_t k(0); k < 3; k++) {
//		n1[k] = c[k] - a[k];
//		d1 -= n1[k] * (a[k] + 0.5 * n1[k]);
//	}
//
//	// swapping planes if necessary
//	if (fabs(n0[0]) < std::numeric_limits<double>::min()) {
//		for (size_t k(0); k < 3; k++)
//			BASELIB::swap(n0[k], n1[k]);
//		BASELIB::swap(d0, d1);
//	}
//
//	if (fabs(n0[0]) < std::numeric_limits<double>::min()) {
//		std::cerr << "points are collinear" << std::endl;
//		std::cout << "E_1: point a: " << a[0] << ", " << a[1] << "," << a[2] <<
//			", point b: " << b[0] << ", " << b[1] << "," << b[2] << ", plane normal: " << n0[0] << ", " << n0[1] << "," << n0[2] << std::endl;
//		std::cout << "E_2: plane normal: " << n1[0] << ", " << n1[1] << "," << n1[2] << std::endl;
//		return;
//	}
//
//	// compute intersection line
//	// g(lambda) = p + lambda * direction
//	// compute the direction vector
//	double direction[3];
//	crossProd(n0, n1, direction);
//
////	std::cout << "direction of intersection line: " << direction[0] << "," << direction[1] << "," << direction[2] << std::endl;
//
//	// compute point in this plane as origin of the intersection line
//	double pnt[3];
//	pnt[1] = (d1 * n0[0] - d0 * n1[0]) / (n0[1] * n1[0] - n1[1] * n0[0]);
//	pnt[0] = -((pnt[1] * n0[1] + d0) / n0[0]);
//	pnt[2] = 0.0; // choose
//
////	std::cout << "point of intersection line: " << pnt[0] << "," << pnt[1] << "," << pnt[2] << std::endl;
//
//	// intersection point between g(lambda) and E_3
//	// describe plane E_3 in normal form (between b and c)
//	double n2[3], d2(0.0);
//	for (size_t k(0); k < 3; k++) {
//		n2[k] = c[k] - b[k];
//		d2 -= n2[k] * (b[k] + 0.5 * n2[k]);
//	}
//
////	std::cout << "E_3: plane normal: " << n2[0] << ", " << n2[1] << "," << n2[2] << std::endl;
//
//	double scpr_direction_n2(scpr(direction, n2, 3));
//	// is intersection line parallel to plane E_3?
//	if (scpr_direction_n2 < std::numeric_limits<double>::min()) {
////		std::cout << "intersection line parallel to E_3" << std::endl;
//		// test if intersection line is in E_3
//		if (scpr(pnt, n2, 3) + d2 < std::numeric_limits<double>::min()) {
//			// middle point of circumscribed sphere
//			for (size_t k(0); k < 3; k++)
//				middle_pnt[k] = pnt[k];
//		} else {
////			std::cerr << "scpr([" << pnt[0] << "," << pnt[1] << "," << pnt[2] << "], [" << n2[0] << "," << n2[1] << "," << n2[2] << ") + " << d2 << ": " << scpr(pnt, n2, 3) + d2 << " "
////					<< "no intersection point between intersection line and plane"
////					<< std::endl;
//			return;
//		}
//	} else {
//		double lambda_s(-(scpr(pnt, n2, 3) + d2) / scpr_direction_n2);
//
//		// middle point of circumscribed sphere
//		for (size_t k(0); k < 3; k++)
//			middle_pnt[k] = pnt[k] + lambda_s * direction[k];
//	}
//
//	// radius of circumscribed sphere
//	sqr_radius = sqrDist(middle_pnt, a);
//}


//bool isEdgeIllegal (const std::vector<GEOLIB::Point*> &pnts, size_t r, LinkedTriangle* triangle)
//{
//	double middle_pnt[3], sqr_radius;
//	getCircumscribedSphereOfTriangle(pnts[(*triangle)[0]]->getData(),
//			pnts[(*triangle)[1]]->getData(), pnts[(*triangle)[2]]->getData(), middle_pnt, sqr_radius);
//
//	if (sqrDist (middle_pnt, pnts[r]->getData()) < sqr_radius) return true;
//
//	return false;
//}

//void legalizeEdge(const std::vector<GEOLIB::Point*> &pnts, size_t r, size_t i,
//		size_t j, LinkedTriangle *triangle)
//{
//	if (isEdgeIllegal (pnts, r, triangle)) {
////		std::cout << "edge " << i << "," << j << " is illegal - flipping " << std::endl;
//		// flip triangle
//		LinkedTriangle *flip_neighbor(NULL);
//		// the neighbor triangles
//		LinkedTriangle *triangle_rnj(NULL), *triangle_lri(NULL);
//		LinkedTriangle *triangle_jpk(NULL), *triangle_ikq(NULL);
//		// point from triangle ijk
//		size_t k;
//
//		// search index k and triangles flip_neighbor and triangle_rnj, triangle_lri
//		if (r == (*triangle)[0]) {
//			// get the triangle, that is in opposite of r
//			flip_neighbor = triangle->getNeighborTriangle(0);
//			if (flip_neighbor == NULL) return;
//			// get node index k from neighbor triangle
//			if ((*flip_neighbor)[0] != i && (*flip_neighbor)[0] != j) k = (*flip_neighbor)[0];
//			else {
//				if ((*flip_neighbor)[1] != i && (*flip_neighbor)[1] != j) k = (*flip_neighbor)[1];
//				else k = (*flip_neighbor)[2];
//			}
//			// get the neighbors of the triangle irj
//			if (i == (*triangle)[1]) {
//				triangle_rnj = triangle->getNeighborTriangle(1);
//				triangle_lri = triangle->getNeighborTriangle(2);
//			} else {
//				triangle_rnj = triangle->getNeighborTriangle(2);
//				triangle_lri = triangle->getNeighborTriangle(1);
//			}
//		}
//		if (r == (*triangle)[1]) {
//			// get the triangle, that is in opposite of r
//			flip_neighbor = triangle->getNeighborTriangle(1);
//			if (flip_neighbor == NULL) return;
//			// get node index k from neighbor triangle
//			if ((*flip_neighbor)[0] != i && (*flip_neighbor)[0] != j)
//				k = (*flip_neighbor)[0];
//			else {
//				if ((*flip_neighbor)[1] != i && (*flip_neighbor)[1] != j)
//					k = (*flip_neighbor)[1];
//				else
//					k = (*flip_neighbor)[2];
//			}
//			// get the neighbors of the triangle irj
//			if (i == (*triangle)[0]) {
//				triangle_rnj = triangle->getNeighborTriangle(0);
//				triangle_lri = triangle->getNeighborTriangle(2);
//			} else {
//				triangle_rnj = triangle->getNeighborTriangle(2);
//				triangle_lri = triangle->getNeighborTriangle(0);
//			}
//		}
//		if (r == (*triangle)[2]) {
//			// get the triangle, that is in opposite of r
//			flip_neighbor = triangle->getNeighborTriangle(2);
//			if (flip_neighbor == NULL) return;
//			// get node index k from neighbor triangle
//			if ((*flip_neighbor)[0] != i && (*flip_neighbor)[0] != j)
//				k = (*flip_neighbor)[0];
//			else {
//				if ((*flip_neighbor)[1] != i && (*flip_neighbor)[1] != j)
//					k = (*flip_neighbor)[1];
//				else
//					k = (*flip_neighbor)[2];
//			}
//			// get the neighbors of the triangle irj
//			if (i == (*triangle)[0]) {
//				triangle_rnj = triangle->getNeighborTriangle(0);
//				triangle_lri = triangle->getNeighborTriangle(1);
//			} else {
//				triangle_rnj = triangle->getNeighborTriangle(1);
//				triangle_lri = triangle->getNeighborTriangle(0);
//			}
//		}
//
//		// get the neighbored triangles of flip_neighbor (opposite of i, j) of the triangle ijk
//		if (i == (*flip_neighbor)[0]) {
//			triangle_jpk = flip_neighbor->getNeighborTriangle(0);
//			if (j == (*flip_neighbor)[1])
//				triangle_ikq = flip_neighbor->getNeighborTriangle(1);
//			else
//				triangle_ikq = flip_neighbor->getNeighborTriangle(2);
//		} else {
//			if (i == (*flip_neighbor)[1]) {
//				triangle_jpk = flip_neighbor->getNeighborTriangle(1);
//				if (j == (*flip_neighbor)[0])
//					triangle_ikq = flip_neighbor->getNeighborTriangle(0);
//				else
//					triangle_ikq = flip_neighbor->getNeighborTriangle(2);
//			} else {
//				triangle_jpk = flip_neighbor->getNeighborTriangle(2);
//				if (j == (*flip_neighbor)[0])
//					triangle_ikq = flip_neighbor->getNeighborTriangle(0);
//				else
//					triangle_ikq = flip_neighbor->getNeighborTriangle(1);
//			}
//		}
//
//		// apply changes to the neighbors of the triangles and vice versa
//
//		// the idx of neighbor of node n in the triangle rnj is the triangle rjk
//		if (triangle_rnj != NULL) {
////			std::cout << "*** rnj *** changing neighbors of triangle: " << *triangle_rnj << std::endl;
////			std::cout << "searching for triangle " << *triangle << " in neighbors: ";
////			triangle_rnj->writeNeighbor (std::cout, 0);
////			triangle_rnj->writeNeighbor (std::cout, 1);
////			triangle_rnj->writeNeighbor (std::cout, 2);
////			std::cout << std::flush;
////			std::cout << "replacing neighbor " << triangle_rnj->getIdxOfNeighborTriangle (triangle) << " by ";
//			triangle_rnj->setNeighborTriangle (triangle_rnj->getIdxOfNeighborTriangle (triangle), flip_neighbor);
////			std::cout << *flip_neighbor << std::endl;
//		}
//
//		// the neighbor of q in the triangle ikq is the triangle irk
//		if (triangle_ikq != NULL) {
////			std::cout << "*** ikq *** changing neighbors of triangle: " << *triangle_ikq << std::endl;
////			std::cout << "searching for triangle " << *flip_neighbor << " in neighbors: ";
////			triangle_ikq->writeNeighbor (std::cout, 0);
////			triangle_ikq->writeNeighbor (std::cout, 1);
////			triangle_ikq->writeNeighbor (std::cout, 2);
////			std::cout << std::endl;
////			std::cout << "replacing neighbor " << triangle_ikq->getIdxOfNeighborTriangle (flip_neighbor) << " by ";
//			triangle_ikq->setNeighborTriangle (triangle_ikq->getIdxOfNeighborTriangle(flip_neighbor), triangle);
////			std::cout << *triangle << std::endl;
//		}
//
//		// the neighbor of r in the triangle irk is the triangle ikq
////		std::cout << "*** irk *** changing neighbors of triangle: " << *triangle << std::endl;
////		std::cout << "searching for triangle " << *flip_neighbor << " in neighbors: ";
////		triangle->writeNeighbor (std::cout, 0);
////		triangle->writeNeighbor (std::cout, 1);
////		triangle->writeNeighbor (std::cout, 2);
////		std::cout << std::endl;
////		std::cout << "replacing neighbor " << triangle->getIdxOfNeighborTriangle (flip_neighbor) << " by ";
//		triangle->setNeighborTriangle (triangle->getIdxOfNeighborTriangle(flip_neighbor), triangle_ikq);
////		std::cout << *triangle_ikq << std::endl;
//
//		// the neighbor of k in the triangle rjk is the triangle rnj
////		std::cout << "k rjk rnj ... " << std::flush;
//		flip_neighbor->setNeighborTriangle (flip_neighbor->getIdxOfNeighborTriangle(triangle), triangle_rnj);
////		std::cout << "ok" << std::endl;
//
//		// the neighbor of i in the triangle irk is the triangle rjk
////		std::cout << "i irk rjk ... " << std::flush;
//		triangle->setNeighborTriangle (triangle->getIdxOfNeighborTriangle(triangle_rnj), flip_neighbor);
////		std::cout << "ok" << std::endl;
//
//		// the neighbor of j in the triangle rjk is the triangle irk
////		std::cout << "j rjk irk ... " << std::flush;
//		flip_neighbor->setNeighborTriangle (flip_neighbor->getIdxOfNeighborTriangle (triangle_ikq), triangle);
////		std::cout << "ok" << std::endl;
//
//		// change triangles indexes
////		std::cout << "old triangle is " << (*triangle) << ": " <<  triangle << std::endl;
////		std::cout << "old flip_neighbor is " << (*flip_neighbor) << ": " <<  flip_neighbor << std::endl;
//		// irj -> irk
//		(*triangle)[triangle->getIdxOfPoint (j)] = k;
//		// ijk -> rjk
//		(*flip_neighbor)[flip_neighbor->getIdxOfPoint (i)] = r;
//
//		legalizeEdge(pnts, r, i, k, triangle);
//		legalizeEdge(pnts, r, k, j, triangle);
//	}
////	else
////		std::cout << "edge " << i << "," << j << " is not illegal" << std::endl;
//}


//void delaunayTriangulation (std::vector<GEOLIB::Point*>& pnts,
//		std::list<LinkedTriangle>& triangles, const Vector& plane_normal)
//{
//	// find the rightmost point among the points with largest y-coordinate
//	size_t idx_rightmost(0);
//	double x_max((*(pnts[idx_rightmost]))[0]), y_max((*(pnts[idx_rightmost]))[1]);
//	size_t n_pnts(pnts.size());
//	for (size_t k(1); k < n_pnts; k++) {
//		if (((*pnts[k]))[0] > x_max) {
//			idx_rightmost = k;
//			x_max = (*(pnts[idx_rightmost]))[0];
//			y_max = (*(pnts[idx_rightmost]))[1];
//		} else {
//			if ( fabs((*(pnts[k]))[0] - x_max) < std::numeric_limits<
//					double>::min() && (*(pnts[k]))[1] > y_max) {
//				idx_rightmost = k;
//				x_max = (*(pnts[idx_rightmost]))[0];
//				y_max = (*(pnts[idx_rightmost]))[1];
//			}
//		}
//	}
//
//	// create two points p_{n} and p_{n+1}, such that the triangle
//	// (p_{n}p_{n+1}pnts[0]) contain all points
//	double x_min_global((*(pnts[0]))[0]);
//	for (size_t k(1); k < n_pnts; k++) {
//		if ((*(pnts[k]))[0] < x_min_global) x_min_global = (*(pnts[k]))[0];
//	}
//
//	double y_min (std::numeric_limits<double>::max());
//	y_max = std::numeric_limits<double>::min();
//	for (size_t k(0); k < n_pnts; k++) {
//		if (k != idx_rightmost) {
//			// calculate line y = m * x + n
//			double m (((*(pnts[idx_rightmost]))[1] - (*(pnts[k]))[1])/((*(pnts[idx_rightmost]))[0] - (*(pnts[k]))[0]));
//			double n ((*(pnts[idx_rightmost]))[1] - m * (*(pnts[idx_rightmost]))[0]);
//			// compute intersection with line at x_min_global
//			double y (m * x_min_global + n);
//			if (y < y_min) y_min = y;
//			if (y > y_max) y_max = y;
//		}
//	}
//
//	// create p_{n} and p_{n+1}
//	GEOLIB::Point pm1 (x_min_global, y_min-(y_max-y_min)*0.001, 0.0);
//	GEOLIB::Point pm2 (x_min_global, y_max+(y_max-y_min)*0.001, 0.0);
//	// check orientation of triangle pnts[idx_rightmost] pm1 pm2
//	Vector u (pm1);
//	u -= Vector (*(pnts[idx_rightmost]));
//	Vector v (pm2);
//	v -= Vector (*(pnts[idx_rightmost]));
//	if (plane_normal.Dot (u.Cross (v)) < 0.0) {
//		pnts.push_back (&pm1);
//		pnts.push_back (&pm2);
//	} else {
//		pnts.push_back (&pm2);
//		pnts.push_back (&pm1);
//	}
//	triangles.push_back (LinkedTriangle (pnts, idx_rightmost, n_pnts, n_pnts+1, NULL, NULL, NULL));
//
//	// check if all points are in the triangle
//	for (size_t r(0); r<n_pnts; r++) {
//		if (! (*(triangles.begin())).containsPoint (*(pnts[r]))) {
//			std::cerr << "initial triangle: " << *(triangles.begin()) << std::endl;
//			std::cerr << *(pnts[idx_rightmost]) << " " << *(pnts[n_pnts]) << " " << *(pnts[n_pnts+1]) << std::endl;
//
//			std::cerr << r << "-th point " << *(pnts[r]) << " not in triangle " << *(triangles.begin()) << std::endl;
//			exit (1);
//		}
//	}
//
//	// delaunay algorithm body starts
//	for (size_t r(0); r<n_pnts; r++) {
//
//		// for checking only
//		if (r != idx_rightmost) {
//			std::cout << "********* triangulation in step " << r << ": " << std::endl;
////			std::list<LinkedTriangle>::iterator it_out (triangles.begin());
////			while (it_out != triangles.end ()) {
////				std::cout << "triangle: (" << (*it_out) << ") with neighbors: ";
////				it_out->writeNeighbor(std::cout,0);
////				std::cout << ", ";
////				it_out->writeNeighbor(std::cout, 1);
////				std::cout << ", ";
////				it_out->writeNeighbor(std::cout, 2);
////				std::cout << std::endl;
////				it_out++;
////			}
//			// checking pointers
//			std::list<LinkedTriangle>::iterator it (triangles.begin());
//			while (it != triangles.end()) {
//				// get neighbors of triangle and check if they have a back pointer
//				for (size_t k(0); k<3; k++) {
//					LinkedTriangle *neighbor (it->getNeighborTriangle (k));
//					if (neighbor) {
//						if (neighbor->getIdxOfNeighborTriangle (&(*it)) == 3) {
//							std::cerr << "triangle " << (*it) << " points to " << *neighbor
//								<< ", but there is no back pointer" << std::endl;
//							exit (1);
//						}
//					}
//				}
//				it++;
//			}
//		}
//
//		if (r != idx_rightmost) {
//			// find the triangle containing p_r
//			std::list<LinkedTriangle>::iterator it (triangles.begin()), it0, it1, it2;
//			bool nfound (true);
//			while (it != triangles.end () && nfound) {
//				if (it->containsPoint (*(pnts[r])))
//					nfound = false;
//				else it++;
//			}
//
//			if (nfound) {
//				std::cerr << "triangle for point " << *(pnts[r]) << " not found" << std::endl;
//				exit (1);
//			} else {
//
//			// split triangle and set neighbors
//			size_t i((*it)[0]), j((*it)[1]), k((*it)[2]);
//			// get neighbors of it
//			LinkedTriangle *neighbor_i (it->getNeighborTriangle(0));
//			LinkedTriangle *neighbor_j (it->getNeighborTriangle(1));
//			LinkedTriangle *neighbor_k (it->getNeighborTriangle(2));
//			it0 = triangles.insert (it, LinkedTriangle(pnts, i,j,r, NULL, NULL, neighbor_k));
//			it1 = triangles.insert (it, LinkedTriangle(pnts, j,k,r, NULL, &(*it0), neighbor_i));
//			it2 = triangles.insert (it, LinkedTriangle(pnts, k,i,r, &(*it0), &(*it1), neighbor_j));
//			// set the missing internal neighbor pointer
//			it0->setNeighborTriangle(0, &(*it1));
//			it0->setNeighborTriangle(1, &(*it2));
//			it1->setNeighborTriangle(0, &(*it2));
//			// set the missing external neighbor pointer
//			if (neighbor_i)
//				neighbor_i->setNeighborTriangle (neighbor_i->getIdxOfNeighborTriangle (&(*it)), &(*it1));
//			if (neighbor_j)
//				neighbor_j->setNeighborTriangle (neighbor_j->getIdxOfNeighborTriangle (&(*it)), &(*it2));
//			if (neighbor_k)
//				neighbor_k->setNeighborTriangle (neighbor_k->getIdxOfNeighborTriangle (&(*it)), &(*it0));
//			// legalize edges
//			legalizeEdge (pnts, r, i, j, &(*it0));
//			legalizeEdge (pnts, r, j, k, &(*it1));
//			legalizeEdge (pnts, r, k, i, &(*it2));
//
//			std::cout << "split triangle [" << *it << "] into [" << *it0 << "], [" << *it1 << "], [" << *it2 << "]" << std::endl;
//
//			// erase splitted triangle
//			triangles.erase (it);
//			}
//		}
//		if (r > 7) r = n_pnts;
//	}
//
//	// remove p_{-1} and p_{-2} and all associated edges
//	// find a triangle containing p_r
//	std::list<LinkedTriangle>::iterator it(triangles.begin());
//	while (it != triangles.end()) {
//		if ((*it)[0] == n_pnts || (*it)[0] == n_pnts+1
//				|| (*it)[1] == n_pnts || (*it)[1] == n_pnts+1
//				|| (*it)[2] == n_pnts || (*it)[2] == n_pnts+1)
//			it = triangles.erase (it);
//		else it++;
//	}
//
//	// remove last two points
//	pnts.erase (--(pnts.end()));
//	pnts.erase (--(pnts.end()));
//}

//void delaunayTriangulatePolygon(const GEOLIB::Polyline* ply, std::list<GEOLIB::Triangle> &triangles)
//{
//	if (!ply->isClosed()) return;
//
//	// copy points
//	size_t n_pnts (ply->getSize()-1);
//	std::vector<GEOLIB::Point*> pnts;
//	pnts.reserve (n_pnts + 2); // + 2 since delaunay algorithm inserts two additional points
//	for (size_t k(0); k < n_pnts; k++)
//		pnts.push_back (new GEOLIB::Point (*(ply->getPoint(k))));
//
////	std::cout << pnts.size() << " initial points: " << std::endl;
////	double x_max ( (*(pnts[0]))[0]), x_min ( (*(pnts[0]))[0] );
////	double y_max ( (*(pnts[0]))[1]), y_min ( (*(pnts[0]))[1] );
////	for (size_t k(0); k < n_pnts; k++) {
////		if ((*(pnts[k]))[0] > x_max) x_max = (*(pnts[k]))[0];
////		if ((*(pnts[k]))[0] < x_min) x_min = (*(pnts[k]))[0];
////		if ((*(pnts[k]))[1] > y_max) y_max = (*(pnts[k]))[1];
////		if ((*(pnts[k]))[1] < y_min) y_min = (*(pnts[k]))[1];
////	}
////	double y_m (1.0 / (y_max-y_min)), x_m (1.0 / (x_max-x_min));
////	if (y_m > x_m) y_m = x_m;
////	else y_m = x_m;
////	double y_n (- y_min * y_m), x_n (- x_min * x_m);
////
////	for (size_t k(0); k < n_pnts; k++) {
////		// std::cout << *(pnts[k]) << " -> ";
////		(*(pnts[k]))[0] = (*(pnts[k]))[0] * x_m + x_n;
////		(*(pnts[k]))[1] = (*(pnts[k]))[1] * y_m + y_n;
////		std::cout << *pnts[k] << " ";
////		std::cout << std::endl;
////	}
////	std::cout << std::endl;
//
//	// rotate copied points into x-y-plane
////	rotatePointsToXY(plane_normal, pnts);
////
////	std::cout << "points after rotate: " << std::endl;
////	for (size_t k(0); k < n_pnts; k++)
////		std::cout << k << " " << *(pnts[k]) << std::endl;
////
////	for (size_t k(0); k < n_pnts; k++) std::cout << k << std::endl;
//
//	// calculate supporting plane
//	Vector plane_normal;
//	double d, z_axis[3];
//	z_axis[0] = 0.0; z_axis[1] = 0.0; z_axis[2] = 1.0;
//
//	getNewellPlane(pnts, plane_normal, d);
//	std::vector<GEOLIB::Point*> pnts_1;
//	pnts_1.reserve (n_pnts);
//
//	size_t orientation;
//	if (scpr (plane_normal.getData(), z_axis, 3) > 0.0) {
//		orientation = CCW;
//		// CCW orientation of polygon
//		for (size_t k(0); k < n_pnts; k++)
//			pnts_1.push_back (new GEOLIB::Point (*(ply->getPoint(k))));
//	} else {
//		orientation = CW;
//		// CW orientation of polygon
//		for (size_t k(0); k < n_pnts; k++) {
//			pnts_1.push_back (new GEOLIB::Point (*(ply->getPoint(n_pnts-1 - k))));
//		}
//	}
//
//	std::list<LinkedTriangle> temp_triangles;
//	delaunayTriangulation (pnts, temp_triangles, plane_normal);
//
//	const std::vector<GEOLIB::Point*>& ref_pnts_vec (ply->getPointsVec());
//	std::list<LinkedTriangle>::const_iterator it (temp_triangles.begin());
//	while (it != temp_triangles.end()) {
//		const size_t i0 (ply->getPointID ((*it)[0]));
//		const size_t i1 (ply->getPointID ((*it)[1]));
//		const size_t i2 (ply->getPointID ((*it)[2]));
//		if (orientation == CCW)
//			triangles.push_back (GEOLIB::Triangle (ref_pnts_vec, i0, i1, i2));
//		else {
//			triangles.push_back (GEOLIB::Triangle (ref_pnts_vec, n_pnts-1-i0, n_pnts-1-i1, n_pnts-1-i2));
//			std::cout << "triangle [" << n_pnts-1-i0 << " " << n_pnts-1-i1 << " " <<  n_pnts-1-i2 << "]" << std::endl;
//		}
//		it++;
//	}
//
//	// delete points
//	for (size_t k(0); k < n_pnts; k++) delete pnts[k];
//	for (size_t k(0); k < n_pnts; k++) delete pnts_1[k];
//}

} // end namespace MATHLIB
