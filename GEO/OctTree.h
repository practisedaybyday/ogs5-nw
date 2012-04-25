/*
 * OctTree.h
 *
 *  Created on: Feb 27, 2012
 *      Author: TF
 */

#ifndef OCTTREE_H_
#define OCTTREE_H_

namespace GEOLIB {

template <typename POINT> class OctTree {
public:
	OctTree(POINT const& ll, POINT const& ur, size_t max_points_per_node) :
		_father (NULL), _ll (ll), _ur (ur), _is_leaf (true),
		_max_points_per_node (max_points_per_node)
#ifndef NDEBUG
		, _depth(0)
#endif
	{
		// init childs
		for (size_t k(0); k < 8; k++) {
			_childs[k] = NULL;
		}

		const double dx(_ur[0] - _ll[0]);
		const double dy(_ur[1] - _ll[1]);
		const double dz(_ur[2] - _ll[2]);

		if (dx > dy && dx > dz) {
			_ll[1] -= (dx-dy)/2.0;
			_ur[1] += (dx-dy)/2.0;
			_ll[2] -= (dx-dz)/2.0;
			_ur[2] += (dx-dz)/2.0;
		} else {
			if (dy > dx && dy > dz) {
				_ll[0] -= (dy-dx)/2.0;
				_ur[0] += (dy-dx)/2.0;
				_ll[2] -= (dy-dz)/2.0;
				_ur[2] += (dy-dz)/2.0;
			} else {
				_ll[0] -= (dz-dx)/2.0;
				_ur[0] += (dz-dx)/2.0;
				_ll[1] -= (dz-dy)/2.0;
				_ur[1] += (dz-dy)/2.0;
			}
		}
//#ifndef NDEBUG
//		std::cout << "root of OctTree: bbx: " << _ll << " x " << _ur << std::endl;
//#endif
	}

	virtual ~OctTree()
	{
		for (size_t k(0); k < 8; k++)
			delete _childs[k];
	}

	/**
	 * This method adds the given point to the OctTree. If necessary,
	 * the OctTree will be extended.
	 * @param pnt the point
	 * @return If the point can be inserted the method returns true, else false.
	 */
	bool addPoint (POINT* pnt)
	{
		if ((*pnt)[0] < _ll[0]) return false;
		if ((*pnt)[0] > _ur[0]) return false;
		if ((*pnt)[1] < _ll[1]) return false;
		if ((*pnt)[1] > _ur[1]) return false;
		if ((*pnt)[2] < _ll[2]) return false;
		if ((*pnt)[2] > _ur[2]) return false;

		if (!_is_leaf) {
			for (size_t k(0); k < 8; k++) {
				if (_childs[k]->addPoint (pnt)) {
					return true;
				}
			}
		}

		// check if point is already in OctTree
		bool pnt_in_tree (false);
		for (size_t k(0); k < _pnts.size() && !pnt_in_tree; k++) {
			const double sqr_dist (MathLib::sqrDist( (_pnts[k])->getData(), pnt->getData() ));
			if (sqr_dist < std::numeric_limits<double>::epsilon())
				pnt_in_tree = true;
		}
		if (!pnt_in_tree)
			_pnts.push_back (pnt);
		else
			return false;

		if (_pnts.size () > _max_points_per_node)
			splitNode ();
		return true;
	}

	/**
	 * range query - returns all points inside the range (min[0], max[0]) x (min[1], max[1]) x (min[2], max[2])
	 * @param min
	 * @param max
	 * @param pnts
	 */
	void getPointsInRange(POINT const& min, POINT const& max, std::vector<POINT*> &pnts) const
	{
		if (_ur[0] < min[0]) return;
		if (_ur[1] < min[1]) return;
		if (_ur[2] < min[2]) return;

		if (max[0] < _ll[0]) return;
		if (max[1] < _ll[1]) return;
		if (max[2] < _ll[2]) return;

		if (_is_leaf) {
			typename std::vector<POINT*>::const_iterator it;
			for (it = (_pnts.begin()); it != _pnts.end(); it++) {
				pnts.push_back(*it);
			}
		} else {
			for (size_t k(0); k<8; k++) {
				_childs[k]->getPointsInRange(min, max, pnts);
			}
		}
	}

private:
	enum OctTreeQuadrant {
		NEL = 0, //!< north east lower
		NWL, //!< north west lower
		SWL, //!< south west lower
		SEL, //!< south east lower
		NEU, //!< south west upper
		NWU, //!< south west upper
		SWU, //!< south west upper
		SEU //!< south east upper
	};

	/**
	 * private constructor
	 * @param ll lower left point
	 * @param ur upper right point
	 * @param father father in the tree
	 * @return
	 */
#ifndef	NDEBUG
	OctTree (POINT const& ll, POINT const& ur, OctTree* father, size_t max_points_per_node, size_t depth) :
#else
	OctTree (POINT const& ll, POINT const& ur, OctTree* father, size_t max_points_per_node) :
#endif
		_father (father), _ll (ll), _ur (ur), _is_leaf (true),
		_max_points_per_node (max_points_per_node)
#ifndef NDEBUG
		, _depth(depth)
#endif
	{
		// init childs
		for (size_t k(0); k < 8; k++)
			_childs[k] = NULL;

//#ifndef NDEBUG
//		if (_depth <= 1)
//			std::cout << "[OctTree] depth: " << _depth << ", bbx: " << _ll << " x " << _ur << std::endl;
//#endif
	}


	void splitNode ()
	{
		const double x_mid((_ur[0] + _ll[0]) / 2.0);
		const double y_mid((_ur[1] + _ll[1]) / 2.0);
		const double z_mid((_ur[2] + _ll[2]) / 2.0);
		POINT p0(x_mid, y_mid, _ll[2]), p1(_ur[0], _ur[1], z_mid);

		// create child NEL
#ifndef NDEBUG
		_childs[NEL] = new OctTree<POINT> (p0, p1, this, _max_points_per_node, _depth+1);
#else
		_childs[NEL] = new OctTree<POINT> (p0, p1, this, _max_points_per_node);
#endif

		// create child NWL
		p0[0] = _ll[0];
		p1[0] = x_mid;
#ifndef NDEBUG
		_childs[NWL] = new OctTree<POINT> (p0, p1, this, _max_points_per_node, _depth+1);
#else
		_childs[NWL] = new OctTree<POINT> (p0, p1, this, _max_points_per_node);
#endif

		// create child SWL
		p0[1] = _ll[1];
		p1[1] = y_mid;
#ifndef NDEBUG
		_childs[SWL] = new OctTree<POINT> (_ll, p1, this, _max_points_per_node, _depth+1);
#else
		_childs[SWL] = new OctTree<POINT> (_ll, p1, this, _max_points_per_node);
#endif

		// create child NEU
#ifndef NDEBUG
		_childs[NEU] = new OctTree<POINT> (p1, _ur, this, _max_points_per_node, _depth+1);
#else
		_childs[NEU] = new OctTree<POINT> (p1, _ur, this, _max_points_per_node);
#endif

		// create child SEL
		p0[0] = x_mid;
		p1[0] = _ur[0];
#ifndef NDEBUG
		_childs[SEL] = new OctTree<POINT> (p0, p1, this, _max_points_per_node, _depth+1);
#else
		_childs[SEL] = new OctTree<POINT> (p0, p1, this, _max_points_per_node);
#endif

		// create child NWU
		p0[0] = _ll[0];
		p0[1] = y_mid;
		p0[2] = z_mid;
		p1[0] = x_mid;
		p1[1] = _ur[1];
		p1[2] = _ur[2];
#ifndef NDEBUG
		_childs[NWU] = new OctTree<POINT> (p0, p1, this, _max_points_per_node, _depth+1);
#else
		_childs[NWU] = new OctTree<POINT> (p0, p1, this, _max_points_per_node);
#endif

		// create child SWU
		p0[1] = _ll[1];
		p1[1] = y_mid;
#ifndef NDEBUG
		_childs[SWU] = new OctTree<POINT> (p0, p1, this, _max_points_per_node, _depth+1);
#else
		_childs[SWU] = new OctTree<POINT> (p0, p1, this, _max_points_per_node);
#endif
		// create child SEU
		p0[0] = x_mid;
		p1[0] = _ur[0];
		p1[1] = y_mid;
		p1[2] = _ur[2];
#ifndef NDEBUG
		_childs[SEU] = new OctTree<POINT> (p0, p1, this, _max_points_per_node, _depth+1);
#else
		_childs[SEU] = new OctTree<POINT> (p0, p1, this, _max_points_per_node);
#endif

		// distribute points to sub quadtrees
		const size_t n_pnts(_pnts.size());
		for (size_t j(0); j < n_pnts; j++) {
			bool nfound(true);
			for (size_t k(0); k < 8 && nfound; k++) {
				if (_childs[k]->addPoint(_pnts[j])) {
					nfound = false;
				}
			}
		}
		_pnts.clear();
		_is_leaf = false;
	}

	OctTree<POINT>* _father;
	/**
	 * childs are sorted:
	 *   _childs[0] is north east lower child
	 *   _childs[1] is north west lower child
	 *   _childs[2] is south west lower child
	 *   _childs[3] is south east lower child
	 *   _childs[4] is north east upper child
	 *   _childs[5] is north west upper child
	 *   _childs[6] is south west upper child
	 *   _childs[7] is south east upper child
	 */
	OctTree<POINT>* _childs[8];
	/**
	 * lower left front face point of the cube
	 */
	POINT _ll;
	/**
	 * upper right back face point of the cube
	 */
	POINT _ur;
	std::vector<POINT*> _pnts;
	bool _is_leaf;
	/**
	 * maximum number of points per leaf
	 */
	const size_t _max_points_per_node;
#ifndef NDEBUG
	size_t _depth;
#endif
};

}

#endif /* OCTTREE_H_ */
