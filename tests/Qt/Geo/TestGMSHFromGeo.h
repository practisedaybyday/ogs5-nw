/*
 * TestGMSHFromGeo.cpp
 * 2012/01/19 KR Initial implementation
 */

#ifndef TESTGMSHFROMGEO_H
#define TESTGMSHFROMGEO_H

#include "QtTestBase.h"
#include "ProjectData.h"
#include "GEOObjects.h"
#include "MeshIO/GMSHInterface.h"
#include "XmlIO/XmlGmlInterface.h"
#include "XmlIO/XmlStnInterface.h"

#include <sstream>


/**
 *
 * @brief Tests the conversion of geometry- / station-data to a *.geo-file
 * needed for mesh-generation with GMSH.
 *
 * What it does:
 * * reads geometry + station files
 * * writes a GMSH geo-file
 *
 * Generated geo-file can be compared to with ground truth.
 *
 **/
class TestGMSHFromGeo : public QtTestBase
{
	Q_OBJECT
	
public:
	/// @brief Constructor which passes argc to the base class constructor.
	TestGMSHFromGeo(int argc) : QtTestBase(argc) {};
	
private slots:
	void testGMSHFromGeo()
	{
		ProjectData project;
		GEOLIB::GEOObjects* geo_objects = new GEOLIB::GEOObjects();
		project.setGEOObjects(geo_objects);

		const std::string source_path (SOURCEPATH);
		std::string schemaFile = source_path + "/FileIO/OpenGeoSysGLI.xsd";
		XmlGmlInterface xml_geo(&project, schemaFile);
		xml_geo.readFile(getTestdataInputDir() + "testGMSHFromGEO-BodeCatchment.gml");
		xml_geo.readFile(getTestdataInputDir() + "testGMSHFromGEO-Bode_Rivers.gml");

		schemaFile = source_path + "/FileIO/OpenGeoSysSTN.xsd";
		XmlStnInterface xml_stn(&project, schemaFile);
		xml_stn.readFile(getTestdataInputDir() + "testGMSHFromGEO-BodeGroundWater.stn");
		xml_stn.readFile(getTestdataInputDir() + "testGMSHFromGEO-BodeSelectedBoreholes.stn");

		std::vector<std::string> geo_names;
		std::vector<std::string> stn_names;
		geo_objects->getGeometryNames(geo_names);
		geo_objects->getStationVectorNames(stn_names);
		for (size_t i=0; i<stn_names.size(); i++)
			geo_names.push_back(stn_names[i]);

		// GMSH parameters
		size_t max_number_of_points_in_quadtree_leaf (10);
		double mesh_density_scaling_pnts(0.5);
		double mesh_density_scaling_stations (0.05);
		
		//std::stringstream.setf(ios::fixed);
		//std::stringstream.precision(6);

		std::string file_name("testGMSHFromGeo.geo");
		FileIO::GMSHInterface gmsh_io(file_name);
		std::string result = gmsh_io.writeAllDataToGMSHInput(*geo_objects, geo_names, max_number_of_points_in_quadtree_leaf, mesh_density_scaling_pnts, mesh_density_scaling_stations);
		
		compareToReference(result.c_str(), QString("testGMSHFromGeo_result.geo"));
	}
};

#endif // TESTGMSHFROMGEO_H