/**
 * \file mapGeoToDEM.cpp
 * 2011/12/19 KR Initial implementation
 */

#include "ProjectData.h"
#include "GEOObjects.h"
#include "FileFinder.h"
#include "VtkRaster.h"
#include "XmlIO/XmlGmlInterface.h"

#include <QApplication>

float* img_data(NULL);	// pixel information
double origin_x(0), origin_y(0), cellsize(0); // image origin + pixel size
size_t width(0), height(0); // image dimensions

float getElevation(size_t x, size_t y)
{
	if ((x<origin_x) || (x>origin_x+(width*cellsize)) || (y<origin_y) || (y>origin_y+(height*cellsize)))
		return 0;

	size_t x_index = static_cast<size_t>((x-origin_x)/cellsize);
	size_t y_index = static_cast<size_t>((y-origin_y)/cellsize);

	return img_data[2*(y_index*width+x_index)];
}

int main (int argc, char* argv[])
{
	QApplication app(argc, argv, false);

	if (argc != 3)
	{
		std::cout << "Changes the z-Coordinates of the geometric objects in the geo-file according to the DEM." << std::endl;
		std::cout << std::endl;
		std::cout << "Usage: " << argv[0] << " <geo-file.gml> <DEM-file.asc>" << std::endl;
		return -1;
	}

	std::string geo_name = argv[1];
	std::string dem_name = argv[2];

	if (geo_name.substr(geo_name.length()-4, 4).compare(".gml") != 0)
	{
		std::cout << "Error: Parameter 1 should be a gml-file" << std::endl;
		std::cout << "Usage: " << argv[0] << " <geo-file.gml> <DEM-file.asc>" << std::endl;
		return -1;
	}

	if (geo_name.substr(geo_name.length()-4, 4).compare(".asc") != 0)
	{
		std::cout << "Error: Parameter 2 should be an asc-file" << std::endl;
		std::cout << "Usage: " << argv[0] << " <geo-file.gml> <DEM-file.asc>" << std::endl;
		return -1;
	}


	ProjectData project;
	GEOLIB::GEOObjects* geo_objects = new GEOLIB::GEOObjects();
	project.setGEOObjects(geo_objects);

	FileFinder fileFinder;
	fileFinder.addDirectory(".");
	fileFinder.addDirectory(std::string(SOURCEPATH).append("/FileIO"));

	XmlGmlInterface xml(&project, fileFinder.getPath("OpenGeoSysGLI.xsd"));
	if (xml.readFile(QString::fromStdString(geo_name)))
	{
		img_data = VtkRaster::loadDataFromASC(dem_name, origin_x, origin_y, width, height, cellsize);

		if (img_data != NULL)
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
			xml.writeFile(new_geo_name, names[0]);

			std::cout << "New file \"" << new_geo_name << " successfully written." << std::endl;
			std::cout << std::endl;
		}

		std::cout << "Error: Could not read DEM file..." << std::endl;
		return -1;
	}
	else
	{
		std::cout << "Error: Could not open geometry file..." << std::endl;
		return -1;
	}
	return -1;
}



