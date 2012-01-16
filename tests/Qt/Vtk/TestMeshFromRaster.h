/*
 * testMeshFromRaster.cpp
 * 2012/01/13 KR Initial implementation
 */

#include "Configure.h"
#include "gtest.h"

#include <iostream>

#include "GridAdapter.h"
#include "VtkMeshConverter.h"
#include "VtkGeoImageSource.h"
#include <vtkSmartPointer.h>
#include <vtkImageData.h>
#include <QtGui/QApplication>

/**
 *
 * Functionality for automated testing
 *
 * What it does:
 * * reads asc-rasterfile into QImage
 * * converts QImage to VtkGeoImageSource
 * * Created OGS-Mesh (GridAdapter) from vtkImageData extracted from VtkGeoImageSource
 * * Writes mesh to file (by converted GridAdapter to CFEMesh)
 *
 * Generated mesh file can be compared to with ground truth.
 *
 **/

TEST(Qt_Vtk, MeshFromRaster)
{
	int argc = 1;
	char** argv;
	argv = new char*[1];
	argv[0] = "testMeshFromRaster";
	QApplication a(argc, &argv[0], false);
	//setlocale(LC_NUMERIC,"C");
	QString fileName(SOURCEPATH);
	fileName += "/UTL/VTK/testMeshFromRaster.asc";
	std::cout << fileName.toStdString() << std::endl;
	vtkSmartPointer<VtkGeoImageSource> geo_image = vtkSmartPointer<VtkGeoImageSource>::New();
	geo_image->setImageFilename(fileName);
	vtkSmartPointer<vtkImageData> image = geo_image->GetOutput();

	GridAdapter* grid = VtkMeshConverter::convertImgToMesh(image, geo_image->getOrigin(), geo_image->getSpacing(), MshElemType::TRIANGLE, UseIntensityAs::ELEVATION);

	std::string result_file_name("testMeshFromRaster_testresult.msh");
	std::ofstream out (result_file_name.c_str());
	if (out.is_open())
		FileIO::OGSMeshIO::write (grid->getCFEMesh(), out);
	out.close();

	// a.exec();
}
