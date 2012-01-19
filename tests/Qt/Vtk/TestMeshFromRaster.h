/*
 * testMeshFromRaster.cpp
 * 2012/01/13 KR Initial implementation
 */

#ifndef TESTMESHFROMRASTER_H
#define TESTMESHFROMRASTER_H

#include "QtTestBase.h"

#include "GridAdapter.h"
#include "VtkMeshConverter.h"
#include "VtkGeoImageSource.h"

#include <sstream>

#include <vtkSmartPointer.h>
#include <vtkImageData.h>

/**
 *
 * @brief Tests the conversion of a raster image to a mesh.
 *
 * What it does:
 * * reads asc-rasterfile into QImage
 * * converts QImage to VtkGeoImageSource
 * * Created OGS-Mesh (GridAdapter) from vtkImageData extracted from VtkGeoImageSource
 * * Compares the created mesh to reference file (by converted GridAdapter to CFEMesh)
 *
 * Generated mesh file can be compared to with ground truth.
 *
 **/
class TestMeshFromRaster : public QtTestBase
{
	Q_OBJECT
	
private slots:
	void test()
	{
		QString fileName = getTestdataInputDir();
		fileName += "testMeshFromRaster.asc";
		vtkSmartPointer<VtkGeoImageSource> geo_image = vtkSmartPointer<VtkGeoImageSource>::New();
		geo_image->setImageFilename(fileName);
		vtkSmartPointer<vtkImageData> image = geo_image->GetOutput();
		image->Update();

		GridAdapter* grid = VtkMeshConverter::convertImgToMesh(image, geo_image->getOrigin(),
			geo_image->getSpacing(), MshElemType::TRIANGLE, UseIntensityAs::ELEVATION);
		
		std::stringstream ss;
		FileIO::OGSMeshIO::write (grid->getCFEMesh(), ss);
		
		compareToReference(ss.str().c_str(), QString("testMeshFromRaster_result.msh"));
	}
};

#endif // TESTMESHFROMRASTER_H