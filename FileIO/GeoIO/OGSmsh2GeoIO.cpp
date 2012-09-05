/*
 * OGSmsh2GeoIO.cpp
 *
 *  Created on: Sep 3, 2012
 *      Author: NW
 */

#include "GeoIO/OGSmsh2GeoIO.h"

#include <vector>

// Base
#include "GEOObjects.h"

// Mesh
#include "msh_mesh.h"
#include "MSHEnums.h"

namespace FileIO
{
void OGSmsh2GeoIO::loadMeshAsGeometry (const MeshLib::CFEMesh* mesh, std::string &project_name, GEOLIB::GEOObjects* geo)
{
    assert(mesh != NULL);
    
	const size_t n_pnts (mesh->getNodeVector().size());
	std::vector<GEOLIB::Point*>* pnts (new std::vector<GEOLIB::Point*>);
	for (size_t k(0); k < n_pnts; k++)
	{
		MeshLib::CNode* node = mesh->nod_vector[k];
		const double* node_xyz = node->getData();
		pnts->push_back (new GEOLIB::Point (node_xyz[0], node_xyz[1], node_xyz[2]));
	}

	geo->addPointVec (pnts, project_name);

	std::vector<size_t> const& pnt_id_map (geo->getPointVecObj(project_name)->getIDMap());
	const size_t n_elements (mesh->getElementVector().size());
	GEOLIB::Surface* sfc (new GEOLIB::Surface (*pnts));
	for (size_t k(0); k < n_elements; k++)
	{
	    const MeshLib::CElem* ele = mesh->ele_vector[k];
		if (ele->GetElementType() == MshElemType::TRIANGLE) // read 3 node triangle
		{
			sfc->addTriangle (pnt_id_map[ele->GetNodeIndex(0)], pnt_id_map[ele->GetNodeIndex(1)], pnt_id_map[ele->GetNodeIndex(2)]);
		}
	}

	std::vector<GEOLIB::Surface*>* sfcs (new std::vector<GEOLIB::Surface*>);
	sfcs->push_back(sfc);
	geo->addSurfaceVec (sfcs, project_name);
}
} // end namespace FileIO
