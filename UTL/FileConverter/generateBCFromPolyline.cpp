/*
 * generateBCFromPolyline.cpp
 *
 *  Created on: Mar 8, 2011
 *      Author: TF
 */

// GEO
#include "GEOObjects.h"
#include "PolylineVec.h"
#include "ProjectData.h"

// FileIO
#include "XMLInterface.h"

#include "problem.h"
Problem* aproblem = NULL;

#include <QString>

int main (int argc, char* argv[])
{
	if (argc == 1)
	{
		std::cout << "Usage: " << argv[0] << " gml-file" << std::endl;
		std::cout << "\tgives you the name of all polylines in the file" << std::endl;
		std::cout << "Usage: " << argv[0] << " gml-file polyline-to-convert" << std::endl;
		std::cout << "\tcreates for the given polyline points boundary conditions" <<
		std::endl;
		return -1;
	}
	GEOLIB::GEOObjects* geo_objs (new GEOLIB::GEOObjects);
	std::string schema_name(
	        "/home/fischeth/workspace/OGS-FirstFloor/sources/FileIO/OpenGeoSysGLI.xsd");
	ProjectData* project_data (new ProjectData);
	project_data->setGEOObjects (geo_objs);
	XMLInterface xml(project_data, schema_name);
	std::string fname (argv[1]);
	xml.readGLIFile(QString::fromStdString (fname));

	std::vector<std::string> geo_names;
	geo_objs->getGeometryNames (geo_names);
	if (geo_names.empty ())
	{
		std::cout << "no geometries found" << std::endl;
		return -1;
	}
	const GEOLIB::PolylineVec* ply_vec (geo_objs->getPolylineVecObj(geo_names[0]));
	if (!ply_vec)
	{
		std::cout << "could not found polylines" << std::endl;
		delete project_data;
		return -1;
	}
	const size_t n_ply (ply_vec->size());

	std::vector<size_t> ply_pnt_ids;
	for (size_t k(0); k < n_ply; k++)
	{
		std::string ply_name;
		if (ply_vec->getNameOfElementByID(k, ply_name))
		{
			if (argc == 2)
				std::cout << "polyline " << k << ": " << ply_name << std::endl;
			else if (ply_name.find (argv[2]) != std::string::npos)
			{
				std::cout << "found polyline " << ply_name << std::endl;
				GEOLIB::Polyline const* ply (ply_vec->getElementByName(ply_name));
				const size_t n_ply_pnts (ply->getNumberOfPoints());
				for (size_t j(0); j < n_ply_pnts; j++)
					ply_pnt_ids.push_back (ply->getPointID(j));
			}
		}
	}

	if (argc == 2)
		return 0;

	std::vector<GEOLIB::Point*> const* geo_pnts (geo_objs->getPointVec(geo_names[0]));
	// write gli file and bc file
	std::ofstream gli_out ("TB.gli");
	std::ofstream bc_out ("TB.bc");
	bc_out << "// file generated by " << argv[0] << std::endl;
	if (gli_out && bc_out)
	{
		gli_out << "#POINTS" << std::endl;
		for (size_t k(0); k < ply_pnt_ids.size(); k++)
		{
			gli_out << k << " " << *((*geo_pnts)[ply_pnt_ids[k]]) << " $NAME PLYPNT" <<
			argv[2] << k << std::endl;
			// boundary condition
			bc_out << "#BOUNDARY_CONDITION" << std::endl;
			bc_out << "\t$PCS_TYPE" << std::endl << "\t\tGROUNDWATER_FLOW" << std::endl;
			bc_out << "\t$PRIMARY_VARIABLE" << std::endl << "\t\tHEAD" << std::endl;
			bc_out << "\t$GEO_TYPE" << std::endl << "\t\tPOINT PLYPNT" << argv[2] <<
			k << std::endl;
			bc_out << "\t$DIS_TYPE" << std::endl << "\t\tCONSTANT " <<
			(*((*geo_pnts)[ply_pnt_ids[k]]))[2] << std::endl;
		}
		gli_out << "#STOP" << std::endl;
		bc_out << "#STOP" << std::endl;
		gli_out.close ();
		bc_out.close ();
	}

	delete project_data;
}

