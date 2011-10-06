/**************************************************************************
MSHLib - Object:
Task:
Programing:
08/2005 OK Encapsulated from mshlib
**************************************************************************/

#include "math.h"
// C++
#include <string>
#include <vector>

// MSHLib
#include "msh_lib.h"
#include "msh_elem.h"
#include "mathlib.h"

// FileIO
#include "MeshIO/GMSHInterface.h"

// PCSLib
#include "rf_pcs.h"
#include "gs_project.h"
#include "rf_mmp_new.h"                           //OK

/**************************************************************************
MSHLib-Method:
Task:
Programing:
07/2005 MB Implementation
11/2005 WW/MB element class
03/2006 OK mat_group
**************************************************************************/
void CFEMesh::CreatePriELEFromTri(int no_layer,double layer_thickness)
{
   no_layer = no_layer;
   layer_thickness = layer_thickness;
   int j;
   int hempel, hampel;
   long i;
   int k;
   long size;
   long no_tri_elements = (long)ele_vector.size();
   long no_tri_nodes = (long)nod_vector.size();
   MeshLib::CElem* m_tri_ele = NULL;
   MeshLib::CElem* m_ele = NULL;
   //----------------------------------------------------------------------
   // Create MSH
   MSHDelete("PRIS_from_TRI");
   CFEMesh* m_msh_pri (new CFEMesh(_geo_obj, _geo_name));
   m_msh_pri->pcs_name = "PRIS_from_TRI";
   m_msh_pri->setElementType (MshElemType::PRISM);
   //m_msh_pri->no_msh_layer = no_layer;
   m_msh_pri->setNumberOfMeshLayers (no_layer);
   //----------------------------------------------------------------------
   // Create Prism elements
   size = (no_layer + 1) * no_tri_nodes;
   m_msh_pri->nod_vector.resize(size);
   for(j=0;j<size;j++)
      m_msh_pri->nod_vector[j] = NULL;
   for(j=0;j<no_layer;j++)
   {
      for(i=0;i<no_tri_elements;i++)
      {
         //Elements
         m_tri_ele = ele_vector[i];
         m_ele = new  MeshLib::CElem;
                                                  //OK
         m_ele->SetPatchIndex((int)mmp_vector.size()+j);
         m_ele->SetElementType(MshElemType::PRISM);
         m_ele->nnodes = 6;
         m_ele->nodes_index.resize(m_ele->nnodes);
         //Set indices
         m_ele->nodes_index[0] = m_tri_ele->GetNodeIndex(0)+ j*no_tri_nodes;
         m_ele->nodes_index[1] = m_tri_ele->GetNodeIndex(1) + j*no_tri_nodes;
         m_ele->nodes_index[2] = m_tri_ele->GetNodeIndex(2) + j*no_tri_nodes;
         m_ele->nodes_index[3] = m_ele->GetNodeIndex(0)+ no_tri_nodes;
         m_ele->nodes_index[4] = m_ele->GetNodeIndex(1)+ no_tri_nodes;
         m_ele->nodes_index[5] = m_ele->GetNodeIndex(2)+ no_tri_nodes;
         //Nodes
         hempel = 0;
         hampel = 0;
         for(k=0;k<m_ele->nnodes;k++)
         {
            if(m_msh_pri->nod_vector[m_ele->GetNodeIndex(k)]==NULL)
            {
               m_ele->nodes[k] = new CNode(m_ele->GetNodeIndex(k));
               m_msh_pri->nod_vector[m_ele->GetNodeIndex(k)] = m_ele->nodes[k];
               if(k>2)
               {
                  hempel = 3;
                  hampel = 1;
               }
               double const*const pnt (nod_vector[m_tri_ele->GetNodeIndex(k-hempel)]->getData());
               m_ele->nodes[k]->SetX(pnt[0]);
               m_ele->nodes[k]->SetY(pnt[1]);
               m_ele->nodes[k]->SetZ(pnt[2] - (j+hampel)*layer_thickness);
            }
            else
            {
               m_ele->nodes[k] = m_msh_pri->nod_vector[m_ele->GetNodeIndex(k)];
            }
         }                                        //end for m_ele->nnodes
         m_msh_pri->ele_vector.push_back(m_ele);
      }                                           //end for no_tri_elements
   }                                              //end for no_layers
   //----------------------------------------------------------------------
   if(m_msh_pri->ele_vector.size()>0)
   {
      m_msh_pri->ConstructGrid();
      fem_msh_vector.push_back(m_msh_pri);
   }
   else
      delete m_msh_pri;
}

/**************************************************************************
MSHLib-Method:
Task:
Programing:
04/2005 OK Implementation
11/2005 MB ELE
**************************************************************************/
void CFEMesh::CreateLineELEFromQuad(int m_numberofprismlayers,double m_thicknessofprismlayer,int m_iMATGroup)
{
   m_iMATGroup = m_iMATGroup;
   int j;
   long i, k;
   double x, y, z;
   long size;
   long no_quad_elements = (long)ele_vector.size();
   //long no_quad_nodes = (long)nod_vector.size();
   MeshLib::CElem* m_quad_ele = NULL;
   MeshLib::CElem* m_ele = NULL;

   //----------------------------------------------------------------------
   // Create MSH
   MSHDelete("LINE_from_QUAD");
   CFEMesh* m_msh_line(new CFEMesh(_geo_obj, _geo_name));
   m_msh_line->pcs_name = "LINE_from_QUAD";
   m_msh_line->setElementType (MshElemType::LINE);
   m_msh_line->setNumberOfMeshLayers (m_numberofprismlayers);
   //----------------------------------------------------------------------
   // Create LINE elements
   size = (m_numberofprismlayers+1) * no_quad_elements;
   m_msh_line->nod_vector.resize(size);

   for(j=0;j<size;j++)
   {
      m_msh_line->nod_vector[j]=NULL;
   }

   //for(j=0;j<m_numberofprismlayers;j++){
   //  for(i=0;i<no_quad_elements;i++){
   for(i=0;i<no_quad_elements;i++)                //Elements
   {
      for(j=0;j<m_numberofprismlayers;j++)
      {
         //Elements
         m_quad_ele = ele_vector[i];
         m_ele = new MeshLib::CElem;
         //m_ele->GetPatchIndex() = m_msh_line->mat_group;
         m_ele->SetIndex((i*m_numberofprismlayers) + j);
         m_ele->SetElementType(MshElemType::LINE);
         m_ele->nnodes = 2;
         m_ele->nodes_index.resize(m_ele->nnodes);
         //Set indices
         m_ele->nodes_index[0] = j + i * (m_numberofprismlayers+1);
         m_ele->nodes_index[1] = m_ele->nodes_index[0]+ 1;
         m_msh_line->ele_vector.push_back(m_ele);
         //
         double const* center(m_quad_ele->GetGravityCenter());
         x = center[0];
         y = center[1];
         z = center[2];
         //Nodes
         for(k=0;k<m_ele->nnodes;k++)
         {
                                                  //if new node
            if(m_msh_line->nod_vector[m_ele->GetNodeIndex(k)]==NULL)
            {
               m_ele->nodes[k] = new CNode(m_ele->GetNodeIndex(k));
               m_msh_line->nod_vector[m_ele->GetNodeIndex(k)] = m_ele->nodes[k];
               //Set Coordinates
               m_ele->nodes[k]->SetX(x);
               m_ele->nodes[k]->SetY(y);
               z = z + ((k+j) * m_thicknessofprismlayer);
               m_ele->nodes[k]->SetZ(z);
            }                                     // end if new node
            else
            {
               m_ele->nodes[k] = m_msh_line->nod_vector[m_ele->GetNodeIndex(k)];
            }
         }                                        // end for nnodes
      }                                           // end for no_layers
   }                                              // end for quad elements
   if(m_msh_line->ele_vector.size()>0)
      fem_msh_vector.push_back(m_msh_line);
   else
      delete m_msh_line;
}

/**************************************************************************
MSHLib-Method:
Task:
Programing:
04/2005 OK Implementation
**************************************************************************/
void CFEMesh::CreateHexELEFromQuad(int no_layer,double layer_thickness)
{
   int j;
   int hempel, hampel;
   long i;
   int k;
   long size;
   long no_quad_elements = (long)ele_vector.size();
   long no_quad_nodes = (long)nod_vector.size();
   MeshLib::CElem* m_quad_ele = NULL;
   MeshLib::CElem* m_ele = NULL;

   //----------------------------------------------------------------------
   // Create MSH
   MSHDelete("HEX_from_QUAD");
   CFEMesh* m_msh_hex (new CFEMesh(_geo_obj, _geo_name));
   m_msh_hex->pcs_name = "HEX_from_QUAD";
   m_msh_hex->setElementType (MshElemType::HEXAHEDRON);
   m_msh_hex->setNumberOfMeshLayers (no_layer);
   //----------------------------------------------------------------------
   // Create HEX elements
   size = (no_layer + 1) * no_quad_nodes;
   m_msh_hex->nod_vector.resize(size);

   for(j=0;j<size;j++)
      m_msh_hex->nod_vector[j] = NULL;

   for(j=0;j<no_layer;j++)
   {
      for(i=0;i<no_quad_elements;i++)
      {
         //Elements
         m_quad_ele = ele_vector[i];
         m_ele = new MeshLib::CElem;
         //m_ele->SetPatchIndex(j);
         m_ele->SetIndex((j*no_quad_elements) + i);
         m_ele->SetElementType(MshElemType::HEXAHEDRON);
         m_ele->nnodes = 8;
         m_ele->nodes_index.resize(m_ele->nnodes);
         //Set indices
         m_ele->nodes_index[0] = m_quad_ele->GetNodeIndex(0) + j*no_quad_nodes;
         m_ele->nodes_index[1] = m_quad_ele->GetNodeIndex(1) + j*no_quad_nodes;
         m_ele->nodes_index[2] = m_quad_ele->GetNodeIndex(2) + j*no_quad_nodes;
         m_ele->nodes_index[3] = m_quad_ele->GetNodeIndex(3) + j*no_quad_nodes;
         m_ele->nodes_index[4] = m_ele->GetNodeIndex(0)+ no_quad_nodes;
         m_ele->nodes_index[5] = m_ele->GetNodeIndex(1)+ no_quad_nodes;
         m_ele->nodes_index[6] = m_ele->GetNodeIndex(2)+ no_quad_nodes;
         m_ele->nodes_index[7] = m_ele->GetNodeIndex(3)+ no_quad_nodes;
         //Nodes
         hempel = 0;
         hampel = 0;
         for(k=0;k<m_ele->nnodes;k++)
         {
            if(m_msh_hex->nod_vector[m_ele->GetNodeIndex(k)]==NULL)
            {
               m_ele->nodes[k] = new CNode(m_ele->GetNodeIndex(k));
               m_msh_hex->nod_vector[m_ele->GetNodeIndex(k)] = m_ele->nodes[k];
               if(k>3)
               {
                  hempel = 4;
                  hampel = 1;
               }
               double const*const pnt (nod_vector[m_quad_ele->GetNodeIndex(k-hempel)]->getData());
               m_ele->nodes[k]->SetX(pnt[0]);
               m_ele->nodes[k]->SetY(pnt[1]);
               m_ele->nodes[k]->SetZ(pnt[2] - (j+hampel)*layer_thickness);
            }
            else
            {
               m_ele->nodes[k] = m_msh_hex->nod_vector[m_ele->GetNodeIndex(k)];
            }
         }                                        //end for m_ele->nnodes
         //msh_no_hexs++;
         m_msh_hex->ele_vector.push_back(m_ele);
      }                                           //end for no_quad_elements
   }                                              //end for no_layers
   //----------------------------------------------------------------------
   if(m_msh_hex->ele_vector.size()>0)
      fem_msh_vector.push_back(m_msh_hex);
   else
      delete m_msh_hex;

}


/**************************************************************************
MSHLib-Method:
Task:
Programing:
04/2005 OK Implementation based on MSH2RFI by WW
08/2005 WW Re-implememtation
10/2005 TK proper ordering and closing of gaps
09/2011 TF changed std::string::compare to std::string::find for the new gmsh format
		to avoid dos unix line-ending issues
**************************************************************************/
void GMSH2MSH(const char* filename,CFEMesh* m_msh)
{
   long id;
   long i = 0;
   int NumNodes = 0;
   int NumElements = 0;
   double x, y, z;
   std::string strbuffer;

   //WW  bool quad=false;
   //WW  CRFProcess* m_pcs = NULL;
   MeshLib::CNode* node = NULL;
   MeshLib::CElem* elem = NULL;
   std::ifstream msh_file(filename, std::ios::in);
   getline(msh_file, strbuffer);                  // Node keyword

   // OLD GMSH  FORMAT----------------------------------------------------------------------
   if (strbuffer.compare("$NOD") == 0)
   {
      while (strbuffer.compare("$ENDELM") != 0)
      {
         msh_file >> NumNodes >> std::ws;
         //....................................................................
         // Node data
         for (i = 0; i < NumNodes; i++)
         {
            msh_file >> id >> x >> y >> z >> std::ws;

            node = new MeshLib::CNode(id, x, y, z);
            m_msh->nod_vector.push_back(node);
         }

         getline(msh_file, strbuffer);            // End Node keyword
         //....................................................................
         // Element data
         getline(msh_file, strbuffer);            // Element keyword
         msh_file >> NumElements >> std::ws;
         for (i = 0; i < NumElements; i++)
         {
            elem = new MeshLib::CElem(i);
            elem->Read(msh_file, 2);
            m_msh->ele_vector.push_back(elem);
         }
         getline(msh_file, strbuffer);            // END keyword

         // ordering nodes and closing gaps TK
         std::vector<int> gmsh_id;
         long new_node_id;
         int counter = 0;
         int diff = 0;
         int j = 0;
         for (i = 0; i < (int) m_msh->nod_vector.size(); i++)
         {
            diff = m_msh->nod_vector[i]->GetIndex() - counter;
            if (diff == 0)
            {
               gmsh_id.push_back(i);
               counter++;
            }
            else
            {
               for (j = 0; j < diff; j++)
               {
                  gmsh_id.push_back(i);
                  counter++;
               }
               i--;
            }
         }

         for (i = 0; i < (int) m_msh->ele_vector.size(); i++)
         {
            for (j = 0; j < (int) m_msh->ele_vector[i]->GetVertexNumber(); j++)
            {
               new_node_id = gmsh_id[m_msh->ele_vector[i]->GetNodeIndex(j)
                  + 1];
               //m_msh->ele_vector[i]->nodes[j]->SetIndex(new_node_id);/*global*/
                                                  /*local*/
               m_msh->ele_vector[i]->nodes_index[j] = new_node_id;
            }
         }
         for (i = 0; i < (int) m_msh->nod_vector.size(); i++)
         {
            m_msh->nod_vector[i]->SetIndex(i);
         }
         // END OF: ordering nodes and closing gaps TK

      }                                           /*End while*/
   }
   // END old GMSH Format----------------------------------------------------------------------
   msh_file.close();


   FileIO::GMSHInterface::readGMSHMesh(filename, m_msh);
//   // NEW 2008 GMSH  FORMAT----------------------------------------------------------------------
//   if (strbuffer.find("$MeshFormat") != std::string::npos)
//   {
//      getline(msh_file, strbuffer);               // version-number file-type data-size
//      getline(msh_file, strbuffer);               //$EndMeshFormat
//      getline(msh_file, strbuffer);               //$Nodes Keywords
//
//      while (strbuffer.find("$EndElements") == std::string::npos)
//      {
//         // Node data
//         msh_file >> NumNodes >> std::ws;
//         for (i = 0; i < NumNodes; i++)
//         {
//            msh_file >> id >> x >> y >> z >> std::ws;
//            node = new MeshLib::CNode(id, x, y, z);
//            m_msh->nod_vector.push_back(node);
//         }
//         getline(msh_file, strbuffer);            // End Node keyword $EndNodes
//
//         // Element data
//         getline(msh_file, strbuffer);            // Element keyword $Elements
//         msh_file >> NumElements >> std::ws;      // number-of-elements
//         for (i = 0; i < NumElements; i++)
//         {
//            elem = new MeshLib::CElem(i);
//            elem->Read(msh_file, 7);
//            if (elem->GetElementType() != MshElemType::INVALID)
//               m_msh->ele_vector.push_back(elem);
//         }
//         getline(msh_file, strbuffer);            // END keyword
//
//         // correct indices TF
//         const size_t n_elements(m_msh->ele_vector.size());
//         for (size_t k(0); k < n_elements; k++)
//         {
//            m_msh->ele_vector[k]->SetIndex(k);
//         }
//
//         // ordering nodes and closing gaps TK
//         std::vector<int> gmsh_id;
//         long new_node_id;
//         int counter = 0;
//         int diff = 0;
//         int j = 0;
//         for (i = 0; i < (int) m_msh->nod_vector.size(); i++)
//         {
//            diff = m_msh->nod_vector[i]->GetIndex() - counter;
//            if (diff == 0)
//            {
//               gmsh_id.push_back(i);
//               counter++;
//            }
//            else
//            {
//               for (j = 0; j < diff; j++)
//               {
//                  gmsh_id.push_back(i);
//                  counter++;
//               }
//               i--;
//            }
//         }
//
//         for (i = 0; i < (int) m_msh->ele_vector.size(); i++)
//         {
//            for (j = 0; j < (int) m_msh->ele_vector[i]->GetVertexNumber(); j++)
//            {
//               new_node_id = gmsh_id[m_msh->ele_vector[i]->GetNodeIndex(j)
//                  + 1];
//               //m_msh->ele_vector[i]->nodes[j]->SetIndex(new_node_id);/*global*/
//                                                  /*local*/
//               m_msh->ele_vector[i]->nodes_index[j] = new_node_id;
//            }
//         }
//         for (i = 0; i < (int) m_msh->nod_vector.size(); i++)
//         {
//            m_msh->nod_vector[i]->SetIndex(i);
//         }
//         // END OF: ordering nodes and closing gaps TK
//
//      }                                           /*End while*/
//   }
//   // END New 2008 GMSH Format----------------------------------------------------------------------

   //  m_msh->ConstructGrid(); // TF


}


/**************************************************************************
MSHLib-Method:
Task:   Takes Surface, meshs it and saves it with defined file name

        surface_name         = name of surface for meshing
        file_name_const_char = path + name without extension

        Attention: You have to remove/delete files extern.

Programing:
12/2005 TK implementation
**************************************************************************/
void Mesh_Single_Surface(std::string surface_name, const char *file_name_const_char)
{

   //Searching Methods
   int i=0, j=0, k=0;
   std::string Name;
   std::vector<CGLPoint*> surface_points_searchvector;
   CGLPoint *m_point = NULL;

   for (i=0; i<(int)surface_points_searchvector.size();i++)
   {
      delete surface_points_searchvector[i];
   }
   surface_points_searchvector.clear();

   //Get SFC
   for (i=0;i<(int)surface_vector.size();i++)
   {
      if (surface_vector[i]->name.find(surface_name)==0)
      {
         for (j=0;j<(int)surface_vector[i]->polyline_of_surface_vector.size();j++)
         {
            if (j==0)
            {
               for (k=0;k<(int)surface_vector[i]->polyline_of_surface_vector[j]->point_vector.size();k++)
               {
                  m_point = new CGLPoint;
                  m_point->nb_of_ply = i;
                  m_point->id = surface_vector[i]->polyline_of_surface_vector[j]->point_vector[k]->id;
                  m_point->x  = surface_vector[i]->polyline_of_surface_vector[j]->point_vector[k]->x;
                  m_point->y  = surface_vector[i]->polyline_of_surface_vector[j]->point_vector[k]->y;
                  m_point->z  = surface_vector[i]->polyline_of_surface_vector[j]->point_vector[k]->z;
                  surface_points_searchvector.push_back(m_point);
                  if (k==(int)surface_vector[i]->polyline_of_surface_vector[j]->point_vector.size()-1 &&
                     surface_vector[i]->polyline_of_surface_vector[j]->point_vector[k]->id == surface_vector[i]->polyline_of_surface_vector[j]->point_vector[0]->id)
                  {
                     surface_points_searchvector.erase(surface_points_searchvector.begin()+k);
                  }
               }
            }

            else
            {
               if (m_point->id != surface_vector[i]->polyline_of_surface_vector[j]->point_vector[0]->id &&
                  m_point->id == surface_vector[i]->polyline_of_surface_vector[j]->point_vector[(int)surface_vector[i]->polyline_of_surface_vector[j]->point_vector.size()-1]->id)
               {
                  for (k=1;k<(int)surface_vector[i]->polyline_of_surface_vector[j]->point_vector.size();k++)
                  {
                     m_point = new CGLPoint;
                     m_point->nb_of_ply = i;
                     m_point->id = surface_vector[i]->polyline_of_surface_vector[j]->point_vector[(int)surface_vector[i]->polyline_of_surface_vector[j]->point_vector.size()-1-k]->id;
                     m_point->x  = surface_vector[i]->polyline_of_surface_vector[j]->point_vector[(int)surface_vector[i]->polyline_of_surface_vector[j]->point_vector.size()-1-k]->x;
                     m_point->y  = surface_vector[i]->polyline_of_surface_vector[j]->point_vector[(int)surface_vector[i]->polyline_of_surface_vector[j]->point_vector.size()-1-k]->y;
                     m_point->z  = surface_vector[i]->polyline_of_surface_vector[j]->point_vector[(int)surface_vector[i]->polyline_of_surface_vector[j]->point_vector.size()-1-k]->z;
                     surface_points_searchvector.push_back(m_point);
                  }
               }
               else
               {
                  for (k=1;k<(int)surface_vector[i]->polyline_of_surface_vector[j]->point_vector.size();k++)
                  {
                     m_point = new CGLPoint;
                     m_point->nb_of_ply = i;
                     m_point->id = surface_vector[i]->polyline_of_surface_vector[j]->point_vector[k]->id;
                     m_point->x  = surface_vector[i]->polyline_of_surface_vector[j]->point_vector[k]->x;
                     m_point->y  = surface_vector[i]->polyline_of_surface_vector[j]->point_vector[k]->y;
                     m_point->z  = surface_vector[i]->polyline_of_surface_vector[j]->point_vector[k]->z;
                     surface_points_searchvector.push_back(m_point);
                  }
               }
            }
         }

         if (surface_points_searchvector[0]->id == surface_points_searchvector[surface_points_searchvector.size()-1]->id)
            surface_points_searchvector.erase(surface_points_searchvector.begin());

         //Write GMSH_GEO_FILE of marked Surface and mesh it
         std::string m_strFileNameGEO = file_name_const_char;
         m_strFileNameGEO = m_strFileNameGEO + ".geo";
         file_name_const_char = m_strFileNameGEO.data();
         FILE *geo_file=NULL;
         geo_file = fopen(file_name_const_char, "w+t");
         long id=0;
         double density = 1e100;
         double topologic_distance;

         for (k=0;k<(int)surface_points_searchvector.size();k++)
         {
            if (k==0)
            {
               density =   EuklVek3dDistCoor ( surface_points_searchvector[k]->x,
                  surface_points_searchvector[k]->y,
                  surface_points_searchvector[k]->z,
                  surface_points_searchvector[k+1]->x,
                  surface_points_searchvector[k+1]->y,
                  surface_points_searchvector[k+1]->z);
               topologic_distance = EuklVek3dDistCoor ( surface_points_searchvector[k]->x,
                  surface_points_searchvector[k]->y,
                  surface_points_searchvector[k]->z,
                  surface_points_searchvector[surface_points_searchvector.size()-1]->x,
                  surface_points_searchvector[surface_points_searchvector.size()-1]->y,
                  surface_points_searchvector[surface_points_searchvector.size()-1]->z);
               if (topologic_distance < density) density = topologic_distance;
            }

            if (k>0 && k<(int)surface_points_searchvector.size()-1)
            {
               density =   EuklVek3dDistCoor ( surface_points_searchvector[k]->x,
                  surface_points_searchvector[k]->y,
                  surface_points_searchvector[k]->z,
                  surface_points_searchvector[k-1]->x,
                  surface_points_searchvector[k-1]->y,
                  surface_points_searchvector[k-1]->z);
               topologic_distance = EuklVek3dDistCoor ( surface_points_searchvector[k]->x,
                  surface_points_searchvector[k]->y,
                  surface_points_searchvector[k]->z,
                  surface_points_searchvector[k+1]->x,
                  surface_points_searchvector[k+1]->y,
                  surface_points_searchvector[k+1]->z);
               if (topologic_distance < density) density = topologic_distance;
            }

            if (k == (int)surface_points_searchvector.size()-1)
            {
               density =   EuklVek3dDistCoor ( surface_points_searchvector[k]->x,
                  surface_points_searchvector[k]->y,
                  surface_points_searchvector[k]->z,
                  surface_points_searchvector[k-1]->x,
                  surface_points_searchvector[k-1]->y,
                  surface_points_searchvector[k-1]->z);
               topologic_distance = EuklVek3dDistCoor ( surface_points_searchvector[k]->x,
                  surface_points_searchvector[k]->y,
                  surface_points_searchvector[k]->z,
                  surface_points_searchvector[0]->x,
                  surface_points_searchvector[0]->y,
                  surface_points_searchvector[0]->z);
               if (topologic_distance < density) density = topologic_distance;
            }

            id++;
            fprintf(geo_file,"%s","Point(");
            fprintf(geo_file,"%li",id);
            fprintf(geo_file,"%s",") = {");
            fprintf(geo_file,"%20.14f",surface_points_searchvector[k]->x);
            fprintf(geo_file,"%s",", ");
            fprintf(geo_file,"%20.14f",surface_points_searchvector[k]->y);
            fprintf(geo_file,"%s",", ");
            fprintf(geo_file,"%20.14f",surface_points_searchvector[k]->z);
            fprintf(geo_file,"%s",", ");
            fprintf(geo_file,"%g",density/4);
            fprintf(geo_file,"%s\n","};");

         }
         id=0;
         for (k=0;k<(int)surface_points_searchvector.size()-1;k++)
         {
            id++;
            fprintf(geo_file,"%s","Line(");
            fprintf(geo_file,"%li",id);
            fprintf(geo_file,"%s",") = {");
            fprintf(geo_file,"%d",k+1);
            fprintf(geo_file,"%s",", ");
            fprintf(geo_file,"%d",k+2);
            fprintf(geo_file,"%s\n","};");

         }
         id++;
         fprintf(geo_file,"%s","Line(");
         fprintf(geo_file,"%li",id);
         fprintf(geo_file,"%s",") = {");
         fprintf(geo_file,"%ld",id);
         fprintf(geo_file,"%s",", 1");
         fprintf(geo_file,"%s\n","};");

         fprintf(geo_file,"%s","Line Loop(");
         fprintf(geo_file,"%li",id+1);
         fprintf(geo_file,"%s",") = {");
         for (k=0;k<(int)surface_points_searchvector.size();k++)
         {
            fprintf(geo_file,"%i",k+1);
            if (k<(int)surface_points_searchvector.size()-1)
               fprintf(geo_file,"%s",", ");
         }
         fprintf(geo_file,"%s\n","};");

         fprintf(geo_file,"%s","Plane Surface(");
         fprintf(geo_file,"%li",id+2);
         fprintf(geo_file,"%s",") = {");
         fprintf(geo_file,"%ld",id+1);
         fprintf(geo_file,"%s\n","};");

         fprintf(geo_file,"%s","Physical Surface(");
         fprintf(geo_file,"%li",id+3);
         fprintf(geo_file,"%s",") = {");
         fprintf(geo_file,"%ld",id+2);
         fprintf(geo_file,"%s\n","};");

         fclose(geo_file);

         std::string m_strExecuteGEO = "gmsh " + m_strFileNameGEO +" -2";

         // PCH & TK: Workaround for the old problem.

         //system(m_strExecute);
         //remove(file_name_const_char);
         //END Meshing
         break;
      }
   }
}

/**************************************************************************
GeoSys-Method:
Task:
Programing:
11/2005 MB
**************************************************************************/
void CFEMesh::SetMSHPart(std::vector<long>&elements_active, long StrangNumber)
{
   int j;
   int k;
   long i;
   long size;
   MeshLib::CElem* m_ele = NULL;
   MeshLib::CFEMesh* m_msh_strang = NULL;
   MeshLib::CNode* m_nod = NULL;
   bool found = false;
   StrangNumber= StrangNumber;
   size = (long)elements_active.size();

   m_msh_strang = FEMGet("MSH_Strang");

   // Create MSH
   if (!m_msh_strang)
   {
      m_msh_strang = new CFEMesh(_geo_obj, _geo_name);
      m_msh_strang->pcs_name = "MSH_Strang";
      m_msh_strang->setElementType (MshElemType::LINE);
      m_msh_strang->setNumberOfMeshLayers (getNumberOfMeshLayers());
      //Resize
      m_msh_strang->ele_vector.resize(size);
      m_msh_strang->Eqs2Global_NodeIndex.resize(size+1);

      fem_msh_vector.push_back(m_msh_strang);
   }

   m_msh_strang->nod_vector.resize(0);

   //Elements
   for(i=0;i<(long)elements_active.size();i++)
   {
      //Elements
      m_ele = ele_vector[elements_active[i]];
      m_msh_strang->ele_vector[i] = m_ele;
      //Nodes
      for(k=0;k<m_ele->nnodes;k++)
      {
         m_nod = m_ele->GetNode(k);
         found = false;
         for(j=0; j<(int)m_msh_strang->nod_vector.size(); j++)
         {
            if(*m_msh_strang->nod_vector[j]==*m_nod)
               found = true;
         }
         if(!found)
            m_msh_strang->nod_vector.push_back(m_nod);
      }
      // cout<< "Element"  << i << "  " << m_ele->GetNodeIndex(0)<< "  "  << m_ele->GetNodeIndex(1) <<  endl;
   }

   for(i=0; i<(long)m_msh_strang->nod_vector.size(); i++)
   {
      m_msh_strang->nod_vector[i]->SetEquationIndex(i);
                                                  //+ test;
      m_msh_strang->Eqs2Global_NodeIndex[i] = m_msh_strang->nod_vector[i]->GetIndex();
      //cout<< " " << i << "  " << m_msh_strang->Eqs2Global_NodeIndex[i] <<  endl;
   }

}

/**************************************************************************
MSHLib-Method:
Task:   CreateLineElementsFromMarkedEdges
Programing:
05/2007 NW implementation
**************************************************************************/
void CFEMesh::CreateLineElementsFromMarkedEdges(CFEMesh*m_msh_ply, std::vector<long> &ele_vector_at_ply)
{
   MeshLib::CElem* m_ele = NULL;
   MeshLib::CEdge* m_edg = NULL;
   Math_Group::vec<CEdge*>ele_edges_vector(15);
   Math_Group::vec<CNode*>edge_nodes(3);
   long no_elements;

   // Create line elements
   std::vector<CEdge*> vct_used_edge;
   for (int i=0;i<(long)ele_vector_at_ply.size();i++)
   {
      m_ele = ele_vector[ele_vector_at_ply[i]];
      m_ele->GetEdges(ele_edges_vector);

      for(int j=0;j<(int)m_ele->GetEdgesNumber();j++)
      {
         m_edg = ele_edges_vector[j];
         if(!m_edg->GetMark())
         {
            continue;
         }
         bool done = false;
         for (int k=0;k<(int)vct_used_edge.size();k++)
         {
            if (vct_used_edge[k] == m_edg)
            {
               done = true;
               break;
            }
         }
         if (done)
         {
            continue;
         }
         else
         {
            vct_used_edge.push_back(m_edg);
         }

         if(m_msh_ply)
         {
            no_elements = (long)m_msh_ply->ele_vector.size();
         }
         else
         {
            no_elements = (long)ele_vector.size();
         }
         CElem* new_ele = new CElem(no_elements);
         new_ele->setElementProperties(MshElemType::LINE);
         m_edg->GetNodes(edge_nodes);
         new_ele->nodes_index[0] = edge_nodes[0]->GetIndex();
         new_ele->nodes_index[1] = edge_nodes[1]->GetIndex();
         new_ele->SetPatchIndex(m_ele->GetPatchIndex()+1);

         if(m_msh_ply)
         {
            m_msh_ply->ele_vector.push_back(new_ele);
         }
         else
         {
            ele_vector.push_back(new_ele);
         }
      }
   }
}
