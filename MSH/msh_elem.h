/**************************************************************************
MSHLib - Object:
Task:
Programing:
08/2005 WW/OK Encapsulation from rf_ele_msh
last modified
**************************************************************************/
#ifndef msh_elem_INC
#define msh_elem_INC

#include <string>
#include <iostream>

// MSHLib
#include "MSHEnums.h"
#include "msh_edge.h"

// PCSLib
namespace process{class CRFProcessDeformation;}
namespace Math_Group{class Matrix;}

class CRFProcess;
//------------------------------------------------------------------------
namespace Mesh_Group
{
   // Process
   using process::CRFProcessDeformation;
   using ::CRFProcess;
   using Math_Group::Matrix;
   //------------------------------------------------------------------------
   // Class definition
   class CElem: public CCore
   {
      public:
         // Methods
         CElem();
         CElem(size_t Index);
                                                  // For Faces: Face, local face index
         CElem(size_t Index, CElem* onwer, int Face);
         CElem(size_t Index, CElem* m_ele_parent);//WWOK

         /**
          * copy constructor
          */
         CElem (CElem const &elem);

         ~CElem();

         //------------------------------------------------------------------

         /**
          * \brief Sets the default properties for the given element type.
          * \param t The element type of this element.
          * \param isFace Signals if the element is initialised as a face of an owner element, some properties of the element might change in this case.
          */
         void setElementProperties(MshElemType::type t, bool isFace = false);

         /**
          * Method (re)sets the patch index. Patch index is used to assign a
          * material to the element.
          * @param pidx an identifier for material
          */
         void setPatchIndex (size_t pidx) { patch_index = pidx; }

         /// Depricated method kept for backward compatibility. Use setElementProperties(MshElemType::type t) instead.
         void Config(MshElemType::type t) { setElementProperties(t); };

         // Geometry
         int GetDimension() const
         {
            return ele_dim;
         }
         double* GetGravityCenter()
         {
            return gravity_center;
         }
         int GetPatchIndex() const
         {
            return patch_index;
         }                                        //MatGroup
         void SetPatchIndex(int value)
         {
            patch_index = value;
         }
         void ComputeVolume();
         void SetFluxArea(double fluxarea)
         {
            area = fluxarea;
         }                                        //CMCD for <3D elements with varying area
         double GetFluxArea()
         {
            return area;
         }                                        //CMCD for <3D elements with varying area

         double calcVolume () const;

         double GetVolume() const
         {
            return volume;
         }
         void SetVolume(double Vol)
         {
            volume = Vol;
         }
         // This will be activated after m_tim->CheckCourant() is ready to work
         // kg44 21042010 activated
         void SetCourant(double Cour)             //CMCD
         {
            courant = Cour;
         }
         double GetCourant()                      //CMCD
         {
            return courant;
         }
         void SetNeumann(double Neum)             //CMCD
         {
            neumann = Neum;
         }
         double GetNeumann()                      //CMCD
         {
            return neumann;
         }
         double GetRepLength()                    //CMCD
         {
            return representative_length;
         }
         //------------------------------------------------------------------
         // ID
         MshElemType::type GetElementType() const { return geo_type; };
         void SetElementType( MshElemType::type type) { geo_type = type; };

         void MarkingAll(bool makop);
         std::string GetName() const;
         //------------------------------------------------------------------
         // Nodes
         vec<long> nodes_index;
         void GetNodeIndeces(vec<long>& node_index) const
         {
            for (int i = 0; i < (int) nodes_index.Size(); i++)
               node_index[i] = nodes_index[i];
         }

         void getNodeIndices(std::vector<size_t>& node_indices)
         {
            for (size_t i = 0; i < nodes_index.Size(); i++)
               node_indices.push_back (nodes_index[i]);
         }

         /**
          * const access to the vector nodes_index
          * @return a const reference to the vector
          */
         const vec<long>& GetNodeIndeces () const { return nodes_index; }

         long GetNodeIndex(int index) const
         {
            return nodes_index[index];
         }
         void SetNodeIndex(int index, long g_index)
         {
            nodes_index[index] = g_index;
         }

         void GetNodes(vec<CNode*>& ele_nodes)
         {
            for (size_t i=0; i<nodes.Size(); i++)
               ele_nodes[i] = nodes[i];
         }

         void GetNodes(std::vector<CNode*>& nodesVec)
         {
            for (size_t i=0; i<nodes.Size(); i++)
               nodesVec.push_back(nodes[i]);
         }

         CNode* GetNode(int index)
         {
            return nodes[index];
         }

         CNode const * GetNode(int index) const
		 {
			return nodes[index];
		 }

         void SetNodes(vec<CNode*>& ele_nodes, bool ReSize = false);
         int GetNodesNumber_H() const
         {
            return nnodesHQ;
         }
         size_t GetNodesNumber(bool quad) const
         {
            if (quad)
               return (size_t)nnodesHQ;
            else
               return (size_t)nnodes;
         }
         int GetVertexNumber() const
         {
            return nnodes;
         }
         void SetNodesNumber(int ivalue)
         {
            nnodes = ivalue;
         }                                        //OK
         CElem* GetOwner()
         {
            return owner;
         }                                        //YD

         // Initialize topological properties
         void InitializeMembers();
         //------------------------------------------------------------------
         // Edges
         void GetEdges(vec<CEdge*>& ele_edges)
         {
            for (int i = 0; i < nedges; i++)
               ele_edges[i] = edges[i];
         }
         CEdge* GetEdge(int index)
         {
            return edges[index];
         }
         void SetEdges(vec<CEdge*>& ele_edges)
         {
            for (int i = 0; i < nedges; i++)
               edges[i] = ele_edges[i];
         }
         // KR not used int FindFaceEdges(const int LocalFaceIndex, vec<CEdge*>& face_edges);
         void SetEdgesOrientation(vec<int>& ori_edg)
         {
            for (int i = 0; i < nedges; i++)
               edges_orientation[i] = ori_edg[i];
         }
         void GetLocalIndicesOfEdgeNodes(const int Edge, int *EdgeNodes);
         int GetEdgesNumber() const
         {
            return nedges;
         }
         //------------------------------------------------------------------
         // Faces
         int GetFacesNumber() const
         {
            return nfaces;
         }

         void SetFace();
         void SetFace(CElem* onwer, const int Face);
         int GetSurfaceFacesNumber() const
         {
            return no_faces_on_surface;
         }
         int GetLocalFaceIndex() const
         {
            return face_index;
         }
         int GetElementFaceNodes(int Face, int *FacesNode);
         //------------------------------------------------------------------

		 // Neighbors
         void SetNeighbors(vec<CElem*>& ele_neighbors)
         {
            for (int i = 0; i < nfaces; i++)
               neighbors[i] = ele_neighbors[i];
         }
         void SetNeighbor(const int LocalIndex, CElem* ele_neighbor)
         {
            neighbors[LocalIndex] = ele_neighbor;
         }
         void GetNeighbors(vec<CElem*>& ele_neighbors)
         {
            for (int i = 0; i < nfaces; i++)
               ele_neighbors[i] = neighbors[i];
         }
         CElem* GetNeighbor(int index)
         {
            return neighbors[index];
         }

		 //------------------------------------------------------------------
         // Coordinates transform
         void FillTransformMatrix();
         void FillTransformMatrix(int noneed);
         double getTransformTensor(int idx);
         void AllocateMeomoryforAngle()
         {
            if (!angle)
               angle = new double[3];
         }                                        // WW
         double GetAngle(int i) const
         {
            return angle[i];
         }                                        // PCH
         void SetAngle(int i, double value)
         {
            angle[i] = value;
         }                                        // PCH
         //------------------------------------------------------------------
         // I/O
         void Read(std::istream& is = std::cin, int fileType = 0);

         void WriteIndex(std::ostream& os = std::cout) const;
         void WriteIndex_TEC(std::ostream& os = std::cout) const;
         void WriteAll(std::ostream& os = std::cout) const;
         void WriteNeighbors(std::ostream& os = std::cout) const;

         //------------------------------------------------------------------
         // MAT
         Vec mat_vector;                          //OKWW
         int matgroup_view;                       //TK
         //------------------------------------------------------------------
         // Operator
         // virtual void operator = (const CElem& elem);
         //-------------------------------------------------------------------

         int selected;
                                                  //YD
         void FaceNormal(const int index0, const int index1, double*);
         double *normal_vector;                   //WW normal_vector[3]; //OK
         void SetNormalVector();                  //OK

         // Since m_tim->CheckCourant() is deactivated, the following member are
         // put in comment.
         // kg44 21042010 reactivated
         double representative_length;            //For stability calculations
         double courant;
         double neumann;                          // MSH topology

      private:
         // Members
         // ID
         MshElemType::type geo_type;              //KR: See MSHEnums.h -  1 Line, 2 Quad, 3 Hex, 4 Tri, 5 Tet, 6 Pris
         CElem* owner;
         // Geometrical properties
         int ele_dim;                             // Dimension of element
         int nnodes;
         int nnodesHQ;
         //
         vec<CNode*> nodes;
         int nedges;
         vec<CEdge*> edges;
         vec<int> edges_orientation;
         int nfaces;
         int no_faces_on_surface;
         int face_index;                          // Local face index for the instance for face
         double volume;
         double gravity_center[3];
         int grid_adaptation;                     // Flag for grid adapting.
         int patch_index;
         /*
          // Since m_tim->CheckCourant() is deactivated, the following member are
          // put in comment.
          double representative_length;//For stability calculations
          double courant;
          double neumann;	  // MSH topology
          */
         double area;                             //Flux area
         //
         // MSH topology
         Matrix* tranform_tensor;
         vec<CElem*> neighbors;
         //vec<CElem*> sons;
         // double angle[3];	// PCH, angle[0] rotation along y axis
         //	    angle[1] rotation along x' axis
         //		angle[2] translation along z'' axis.
         double *angle;                           // Dymanic allocate memory.  WW
         //WW double MatT[9];

         // -- Methods
         int GetElementFaces1D(int *FaceNode);
         int GetElementFacesTri(int Face, int *FaceNode);
         int GetElementFacesQuad(int Face, int *FaceNode);
         int GetElementFacesHex(int Face, int *FaceNode);
         int GetElementFacesTet(int Face, int *FaceNode);
         int GetElementFacesPri(int Face, int *FaceNode);
         //-- Friends
         friend class CFEMesh;
         // FEM
         friend class FiniteElement::CElement;
         friend class FiniteElement::CFiniteElementStd;
         friend class FiniteElement::CFiniteElementVec;
         friend class FiniteElement::ElementMatrix;
         friend class FiniteElement::ElementMatrix_DM;
         // PCS
         friend class process::CRFProcessDeformation;
         friend class ::CRFProcess;
   };

} // namespace Mesh_Group
#endif
