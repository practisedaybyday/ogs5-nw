/**
 * \file mapGeoToMesh.cpp
 * 2011/12/19 KR Initial implementation
 */

#include "ProjectData.h"
#include "GEOObjects.h"
#include "FileFinder.h"
#include "OGSMeshIO.h"
#include "MshEditor.h"
#include "XmlIO/XmlGmlInterface.h"
#include "XmlIO/XmlStnInterface.h"

#include <QApplication>
#include <QString>

#include <fstream>
#include <iostream>

MeshLib::CFEMesh* mesh;
MeshLib::CFEMesh flat_sfc_mesh;

float getElevation(size_t x, size_t y)
{
	long idx = flat_sfc_mesh.GetNODOnPNT(new GEOLIB::Point(x,y,0));

	MeshLib::CNode* node = mesh->nod_vector[flat_sfc_mesh.nod_vector[idx]->GetIndex()];

	return node->getData()[2];
}

int mapGeometry(const std::string &geo_name)
{
	std::cout << "Mapping " << geo_name << std::endl;
	ProjectData project;
	GEOLIB::GEOObjects* geo_objects = new GEOLIB::GEOObjects();
	project.setGEOObjects(geo_objects);

	FileFinder fileFinder;
	fileFinder.addDirectory(".");
	fileFinder.addDirectory(std::string(SOURCEPATH).append("/FileIO"));

	FileIO::XmlGmlInterface xml(&project, fileFinder.getPath("OpenGeoSysGLI.xsd"));
	if (xml.readFile(QString::fromStdString(geo_name)))
	{
		std::vector<std::string> names;
		geo_objects->getGeometryNames(names);
		std::vector<GEOLIB::Point*> *points = const_cast<std::vector<GEOLIB::Point*>*>(geo_objects->getPointVec(names[0]));
		size_t nPoints (points->size());
		for (size_t j=0; j<nPoints; j++)
		{
			GEOLIB::Point* pnt = (*points)[j];
			(*pnt)[2] = getElevation((*pnt)[0],(*pnt)[1]);
		}

		std::string new_geo_name = geo_name.substr(0, geo_name.length()-4) + "_elevation.gml";
		xml.setNameForExport(names[0]);
		xml.writeToFile(new_geo_name);

		std::cout << "New file \"" << new_geo_name << " successfully written." << std::endl;
		std::cout << std::endl;

		return 1;
	}

	std::cout << "Error: Could not open geometry file..." << std::endl;
	return -1;
}

int mapStations(const std::string &geo_name)
{
	std::cout << "Mapping " << geo_name << std::endl;
	ProjectData project;
	GEOLIB::GEOObjects* geo_objects = new GEOLIB::GEOObjects();
	project.setGEOObjects(geo_objects);

	FileFinder fileFinder;
	fileFinder.addDirectory(".");
	fileFinder.addDirectory(std::string(SOURCEPATH).append("/FileIO"));

	FileIO::XmlStnInterface xml(&project, fileFinder.getPath("OpenGeoSysSTN.xsd"));
	if (xml.readFile(QString::fromStdString(geo_name)))
	{
		bool is_borehole(false);
		std::vector<std::string> names;
		geo_objects->getStationVectorNames(names);
		std::vector<GEOLIB::Point*> *points = const_cast<std::vector<GEOLIB::Point*>*>(geo_objects->getStationVec(names[0]));
		if (static_cast<GEOLIB::StationBorehole*>((*points)[0])->type() == GEOLIB::Station::BOREHOLE)
			is_borehole = true;
		size_t nPoints (points->size());
		for (size_t j=0; j<nPoints; j++)
		{
			GEOLIB::Point* pnt = (*points)[j];

			if (!is_borehole)
				(*pnt)[2] = getElevation((*pnt)[0],(*pnt)[1]);
			else
			{
				double offset(getElevation((*pnt)[0],(*pnt)[1]) - (*pnt)[2]);
				GEOLIB::StationBorehole* borehole = static_cast<GEOLIB::StationBorehole*>(pnt);
				const std::vector<GEOLIB::Point*> layers = borehole->getProfile();
				size_t nLayers = layers.size();
				for (size_t k=0; k<nLayers; k++)
				{
					GEOLIB::Point* layer_pnt = layers[k];
					(*layer_pnt)[2] = (*layer_pnt)[2] + offset;
				}
			}
		}
	
		std::string new_geo_name = geo_name.substr(0, geo_name.length()-4) + "_elevation.stn";
		xml.setNameForExport(names[0]);
		xml.writeToFile(new_geo_name);

		std::cout << "New file \"" << new_geo_name << " successfully written." << std::endl;
		std::cout << std::endl;

		return 1;
	}

	std::cout << "Error: Could not open geometry file..." << std::endl;
	return -1;
}

void constructFlatSurfaceMesh()
{
	//MshEditor::getSurfaceNodes(*mesh);
	MeshLib::CFEMesh* sfc_mesh = MshEditor::getMeshSurface(*mesh);

	size_t nNodes = sfc_mesh->nod_vector.size();
	for (size_t i=0; i<nNodes; i++)
	{
		const double* coords = sfc_mesh->nod_vector[i]->getData();
		flat_sfc_mesh.nod_vector.push_back(new MeshLib::CNode(coords[0],coords[1],0));
		flat_sfc_mesh.nod_vector[i]->SetIndex(sfc_mesh->nod_vector[i]->GetIndex());
	}
	size_t nElems = sfc_mesh->ele_vector.size();
	for (size_t i=0; i<nElems; i++)
	{
		MeshLib::CElem* old_elem = sfc_mesh->ele_vector[i];
		MeshLib::CElem* elem = new MeshLib::CElem(MshElemType::TRIANGLE, old_elem->GetNodeIndex(0), old_elem->GetNodeIndex(1), old_elem->GetNodeIndex(2), 0);
		flat_sfc_mesh.ele_vector.push_back(elem);
	}
	flat_sfc_mesh.ConstructGrid();
}

int main (int argc, char* argv[])
{
	/*
	QApplication app(argc, argv, false);

	if (argc != 3)
	{
		std::cout << "Changes the z-Coordinates of the geometric objects in the gml- or stn-files according to a given Mesh." << std::endl;
		std::cout << std::endl;
		std::cout << "Usage: " << argv[0] << " <geo-file.gml> <DEM-file.asc>" << std::endl;
		return -1;
	}
	
	bool isList(false);
	std::string geo_name  = argv[1];
	std::string mesh_name = argv[2];
	std::string gml_name("");

	if ((geo_name.substr(geo_name.length()-4, 4).compare(".gml") != 0) && 
		(geo_name.substr(geo_name.length()-4, 4).compare(".stn") != 0) &&
		(geo_name.substr(geo_name.length()-4, 4).compare(".lst") != 0)) 
	{
		std::cout << "Error: Parameter 1 should be a gml- or stn-file" << std::endl;
		std::cout << "Usage: " << argv[0] << " <geo-file.gml> <DEM-file.asc>" << std::endl;
		std::cout << std::endl;
		return -1;
	}

	if (mesh_name.substr(mesh_name.length()-4, 4).compare(".msh") != 0)
	{
		std::cout << "Error: Parameter 2 should be a mesh-file" << std::endl;
		std::cout << "Usage: " << argv[0] << " <geo-file.gml> <mesh-file.msh>" << std::endl;
		return -1;
	}

	if (geo_name.substr(geo_name.length()-4, 4).compare(".lst") == 0)
		isList = true;
	*/
	bool isList(false);
	std::string geo_name  = "c:/project/data/ammer/WESSRivers.gml";
	std::string mesh_name = "c:/project/data/ammer/Ammer-Homogen100m-Final.msh";
	std::string gml_name("");

	FileIO::OGSMeshIO mesh_io;
	mesh = mesh_io.loadMeshFromFile(mesh_name);

	if (mesh != NULL)
	{
		constructFlatSurfaceMesh();

		// map list of geometries to the same DEM
		if (isList)
		{
			std::ifstream in(geo_name.c_str());
			while (!in.eof())
			{
				in >> gml_name;

				if (gml_name.substr(gml_name.length()-4, 4).compare(".gml") == 0)
					mapGeometry(gml_name);
				else if (gml_name.substr(gml_name.length()-4, 4).compare(".stn") == 0)
					mapStations(gml_name);
				else
					std::cout << "File extension for " << gml_name << " unknown." << std::endl;
			}
			return 1;
		}
		// map only one geometry / station file
		else
		{
			if (geo_name.substr(geo_name.length()-4, 4).compare(".gml") == 0)
				mapGeometry(geo_name);
			else if (geo_name.substr(geo_name.length()-4, 4).compare(".stn") == 0)
				mapStations(geo_name);
		}
	}
	else
	{
		std::cout << "Error: Could not read mesh file..." << std::endl;
		return -1;
	}
	return -1;
}



