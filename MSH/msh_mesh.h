/**
 * \file msh_mesh.h
 */

/**************************************************************************
 MSHLib - Object:
 Task:
 Programing:
 08/2005 WW/OK Encapsulation from rf_ele_msh
 last modified
 **************************************************************************/
#ifndef msh_mesh_INC
#define msh_mesh_INC

// C++
#include "Configure.h"
#include <string>

/** depreciated includes of GEOLib */
#include "geo_lib.h"
//#include "geo_pnt.h"
#include "geo_sfc.h"
#include "geo_vol.h"

// GEOLIB
#include "GEOObjects.h"
#include "Point.h"
#include "Polyline.h"
#include "Surface.h"

// MSHLib
#include "MSHEnums.h" // KR 2010/11/15
#include "msh_elem.h"
#include "MeshNodesAlongPolyline.h"

class RandomWalk;
class CFluidMomentum;

#ifdef NEW_EQS    //1.11.2007 WW
namespace Math_Group {class SparseTable;}
using Math_Group::SparseTable;
#endif


//------------------------------------------------------------------------
namespace Mesh_Group {

//------------------------------------------------------------------------
// Class definition
class CFEMesh {

public:
	/// Constructor using geometric information.
	CFEMesh(GEOLIB::GEOObjects* geo_obj = NULL, std::string* unique_name = NULL);
	
	/// Copy-Constructor.
	CFEMesh(const CFEMesh* mesh);

	/// Destructor
	~CFEMesh();

	GEOLIB::GEOObjects* getGEOObjects () const { return _geo_obj; }

	std::string* getProjectName () const { return _geo_name; }

	/**
	 * sets the value for element type
	 * @param ele_type
	 */
	void setElementType (MshElemType::type ele_type);

	/**
	 * set the number of mesh layer
	 * @param n_msh_layer
	 */
	void setNumberOfMeshLayers (size_t n_msh_layer);

	/**
	 * returns the number of mesh layers
	 * @return the number of mesh layers
	 */
	size_t getNumberOfMeshLayers () const;

	/**
	 *
	 * @return
	 */
	bool hasCrossSection() const;

	size_t getNumberOfLines () const;
	size_t getNumberOfQuads () const;
	size_t getNumberOfHexs () const;
	size_t getNumberOfTris () const;
	size_t getNumberOfTets () const;
	size_t getNumberOfPrisms () const;
	double getMinEdgeLength () const;
	/**
	 * do not use this method REMOVE CANDIDATE
	 * @param val
	 */
	void setMinEdgeLength (double val) {_min_edge_length = val;} // TF temporary

	ios::pos_type Read(std::ifstream*);

#ifdef USE_TOKENBUF
	int Read(TokenBuf* tokenbuf);
#endif

	void Write(fstream*) const;
	ios::pos_type GMSReadTIN(std::ifstream*);
	//
	void ConstructGrid();
	void GenerateHighOrderNodes();
	//
	void RenumberNodesForGlobalAssembly();
	// For number of nodes
	int GetMaxElementDim() {
		return max_ele_dim;
	}
	void SwitchOnQuadraticNodes(bool quad) {
		useQuadratic = quad;
	}
	bool getOrder() const {
		return useQuadratic;
	}
	bool isAxisymmetry() const {
		return _axisymmetry;
	}
	// Get number of nodes
	long GetNodesNumber(const bool quadr) const//CMCD int to long
	{
		if (quadr)
			return NodesNumber_Quadratic;
		else
			return NodesNumber_Linear;
	}

	long NodesInUsage() const {
		if (useQuadratic)
			return NodesNumber_Quadratic;
		else
			return NodesNumber_Linear;
	}

	void InitialNodesNumber()//WW
	{
		NodesNumber_Quadratic = NodesNumber_Linear = nod_vector.size();
	}

	/**
	 * @{
	 * */

	/**
	 * \callgraph
	 * \callergraph
	 * \brief depreciated method
	 * @param ply CGLPolyline
	 * @param msh_nod_vector
	 */
	void GetNODOnPLY(CGLPolyline* ply, std::vector<long>& msh_nod_vector);
	/**
	 * \brief depreciated method
	 */
	void GetNodesOnArc(CGLPolyline*m_ply, std::vector<long>&msh_nod_vector);

	/**
	 * \brief depreciated method - uses old surface class
	 */
	void GetNODOnPLY_XY(CGLPolyline*m_ply, std::vector<long>&msh_nod_vector);
//	/**
//	 * \brief depreciated method - uses old surface class
//	 */
	void GetNODOnSFC_Vertical(Surface*m_sfc,vector<long>&msh_nod_vector);

	/**
	 * \brief depreciated method
	 */
	void CreateLineELEFromPLY(CGLPolyline*, int, CFEMesh*);
	/**
	 * \brief depreciated method
	 */
	void CreateLineELEFromPLY(CGLPolyline*);
	/**
	 * \brief depreciated method
	 */
	void CreateLayerPolylines(CGLPolyline*); //OK
	/**
	 * \brief depreciated method
	 */
	void GetELEOnPLY(CGLPolyline*, std::vector<long>&); //OK

	// GEO-SFC
	/**
	 * \brief depreciated method
	 */
	void GetNODOnSFC(Surface*, std::vector<long>&);
	/**
	 * \brief depreciated method
	 */
	void GetNODOnSFC_PLY(Surface*, std::vector<long>&);
	/**
	 * \brief depreciated method
	 */
	void GetNODOnSFC_PLY_XY(Surface*m_sfc, std::vector<long>& msh_nod_vector,
			bool givenNodesOnSurface = false); // givenNodeOnSurface by WW
	/**
	 * \brief depreciated method
	 */
	void GetNODOnSFC_PLY_Z(Surface*,vector<long>&); // 02.2009/OK
	/**
	 * \brief depreciated method
	 */
	void GetNODOnSFC_TIN(Surface*, std::vector<long>&);
	/**
	 * \brief depreciated method
	 */
	void GetNodesOnCylindricalSurface(Surface*m_sfc, std::vector<long>& NodesS);
	/**
	 * \brief depreciated method
	 */
	void CreateQuadELEFromSFC(Surface*);
	/**
	 * \brief depreciated method
	 */
	void CopySelectedNodes(std::vector<long>&msh_nod_vector); //TK

	// REMOVE CANDIDATE
	/**
	 * \brief depreciated method
	 */
	void CreateLineELEFromSFC(); //OK
	// GEO-VOL

	/**
	 * GetNODOnPNT searchs the nearest node to the geometric point
	 * */
	long GetNODOnPNT(const GEOLIB::Point* const pnt) const;
	/**
	 * GetNearestELEOnPNT searchs the nearest element (gravity center)
	 * to the geometric point
	 * */
	long GetNearestELEOnPNT(const GEOLIB::Point* const pnt) const;

	/**
	 * GetNODOnPLY search the nearest nodes along the Polyline object
	 * @param ply constant pointer to a constant Polyline object
	 * @param msh_nod_vector the mesh node indices are saved in this vector
	 * */
	void GetNODOnPLY(const GEOLIB::Polyline* const ply,
			std::vector<size_t>& msh_nod_vector);

	/**
	 *
	 * @param ply
	 * @return
	 */
	const MeshNodesAlongPolyline& GetMeshNodesAlongPolyline(const GEOLIB::Polyline* const ply);

	/**
	 *
	 * @param ply
	 * @param points
	 */
	void getPointsForInterpolationAlongPolyline (const GEOLIB::Polyline* const ply, std::vector<double>& points);

	/**
	 * GetNODOnPLY search the nearest nodes to the Polyline
	 * */
	void GetNODOnPLY(const GEOLIB::Polyline* const ply, std::vector<long>& msh_nod_vector);

	/**
	 * \brief get nodes near the circle arc described by the middle point m, the arc start point a
	 * and the arc end point b.
	 *
	 * If the angle is to small (a == b) then all mesh nodes within the annulus defined by
	 * the inner radius \f$ \|(a-m) \| - min\_edge\_length \f$ and the outer radius
	 * \f$\|(a-m) \| + min\_edge\_length \f$ are pushed in msh_nod_vector
	 */
	void GetNodesOnArc(const GEOLIB::Point* a, const GEOLIB::Point* m,
			const GEOLIB::Point* b, std::vector<size_t>& msh_nod_vector) const;

	/**
	 * \brief gives the indices of CElement elements, which have an edge
	 * in common with the polyline.
	 */
	void GetELEOnPLY(const GEOLIB::Polyline*, std::vector<size_t>&);

	/**
	 * \brief gives the indices of nodes, which are contained in the surface
	 */
	void GetNODOnSFC(const GEOLIB::Surface* sfc, std::vector<size_t>& msh_nod_vector) const;

	/** @} */ // close doxygen group

	//....................................................................
	// QUAD->HEX
	void CreateHexELEFromQuad(int, double);
	// QUAD->LINE
	void CreateLineELEFromQuad(int, double, int);
	void SetActiveElements(std::vector<long>&);
	void SetMSHPart(std::vector<long>&, long); //MB
	bool NodeExists(size_t node);
	// NOD-ELE relations
	//OK411 void SetNOD2ELETopology();
	// ELE-NOD relations
	void SetELE2NODTopology();
	// LINE->LINE
	void AppendLineELE();
	// TRI->PRIS
	void CreatePriELEFromTri(int, double);
	// TRI->LINE
	void CreateLineELEFromTri(); //OK
	void CreateLineELEFromTriELE(); //OK


	// Coordinate system
	int GetCoordinateFlag() const { return coordinate_system; };
	void FillTransformMatrix();

	/**
	 * returns the vector storing pointers to all nodes (class CNode) of the mesh
	 * @return
	 */
	const std::vector<Mesh_Group::CNode*>& getNodeVector () const { return nod_vector; }

	// All nodes - should be private!!!
	std::vector<Mesh_Group::CNode*> nod_vector;

	// All edges
	std::vector<Mesh_Group::CEdge*> edge_vector;
	// All surface feces
	std::vector<Mesh_Group::CElem*> face_vector;
	// All surface nomal
	std::vector<double*> face_normal; //YD
	/**
	 * all elements stored in this vector
	 * */
	std::vector<Mesh_Group::CElem*> ele_vector;
	// Nodes in usage
	// To record eqs_index->global node index
	std::vector<long> Eqs2Global_NodeIndex;

	void PrismRefine(const int Layer, const int subdivision); //OK
	//	void EdgeLengthMinMax(); //OK
	// TF the following two methods are not used, at least in the standard config
//	void SetMATGroupFromVOLLayer(CGLVolume*); //OK
//	void SetMATGroupsFromVOLLayer(); //OK
	void ConnectedNodes(bool quadratic) const; //OK
	void ConnectedElements2Node(bool quadratic = false); //WW
//	void CreateSurfaceTINfromTri(Surface*); //OK
//	void CreateLayerSurfaceTINsfromPris(Surface*); //OK
	// MAT
	std::vector<std::string> mat_names_vector; //OK
	void DefineMobileNodes(CRFProcess*); //OK
	void FaceNormal(); // YD
	void SetNODPatchAreas(); //OK4310
	void SetNetworkIntersectionNodes(); //OK4319->PCH
#ifdef NEW_EQS   // 1.11.2007 WW
	// Compute the graph of the sparse matrix related to this mesh. 1.11.2007 WW
	void CreateSparseTable();
	// Get the sparse graph   1.11.2007 WW
	SparseTable *GetSparseTable(bool quad = false)
	const {if(!quad) return sparse_graph; else return sparse_graph_H;}
#endif

	std::string pcs_name;
	std::string geo_name; //MB
	std::string geo_type_name; //OK10_4310

	int max_mmp_groups; //OKCC
	int msh_max_dim;

	RandomWalk* PT; // PCH

	CFluidMomentum* fm_pcs; // by PCH


private:
	// private attributes
	/**
	 * reference to object of class GEOObject, that manages the geometry data
	 */
	GEOLIB::GEOObjects* _geo_obj;
	/**
	 * identifier for geometry
	 */
	std::string* _geo_name;

	MshElemType::type _ele_type;
	size_t _n_msh_layer; //OK
	bool _cross_section;
	size_t _msh_n_lines;
	size_t _msh_n_quads;
	size_t _msh_n_hexs;
	size_t _msh_n_tris;
	size_t _msh_n_tets;
	size_t _msh_n_prisms;

	double _min_edge_length; //TK

	std::vector<Mesh_Group::MeshNodesAlongPolyline> _mesh_nodes_along_polylines;


	// Process friends
	friend class ::CRFProcess;

	size_t NodesNumber_Linear;
	size_t NodesNumber_Quadratic;
	bool useQuadratic;
	bool _axisymmetry;

	// Coordinate indicator
	// 1:  X component only
	// 12: Y component only
	// 13: Z component only
	// 2:  X, Y component
	// 23:  X, Z component
	// 3:  X, Y, Z component
	int coordinate_system;
	bool has_multi_dim_ele;
	int max_ele_dim;
	int map_counter; //21.01.2009 WW
	bool mapping_check; //23.01.2009 WW

	// Sparse graph of this mesh. 1.11.2007 WW
#ifdef NEW_EQS
	SparseTable *sparse_graph;
	SparseTable *sparse_graph_H; // For high order interpolation
#endif
	// LINE
	void CheckMarkedEdgesOnPolyLine(CGLPolyline*m_polyline,
			std::vector<long> &ele_vector_at_ply); //NW
	void CreateLineElementsFromMarkedEdges(CFEMesh*m_msh_ply,
			std::vector<long> &ele_vector_at_ply); //NW
	bool HasSameCoordinatesNode(CNode* nod, long &node_no); //NW
};

} // namespace Mesh_Group

#endif
