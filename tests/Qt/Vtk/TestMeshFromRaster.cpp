/*
 * testMeshFromRaster.cpp
 * 2012/01/13 KR Initial implementation
 */

#include "gtest.h"
#include "TestHelperFunctions.h"

#include "GridAdapter.h"
#include "VtkMeshConverter.h"
#include "VtkGeoImageSource.h"

#include <sstream>

#include <vtkSmartPointer.h>
#include <vtkImageData.h>

TEST(VTK, TestMeshFromRaster)
{
	QString fileName = getTestdataInputDir();
	fileName += "testMeshFromRaster.asc";
	vtkSmartPointer<VtkGeoImageSource> geo_image = vtkSmartPointer<VtkGeoImageSource>::New();
	geo_image->setImageFilename(fileName);
	vtkSmartPointer<vtkImageData> image = geo_image->GetOutput();
	image->Update();
	
	double origin[3];
	geo_image->getOrigin(origin[0], origin[1], origin[2]);

	GridAdapter* grid = VtkMeshConverter::convertImgToMesh(image, origin, geo_image->getSpacing(), MshElemType::TRIANGLE, UseIntensityAs::ELEVATION);
		
	// Correct number of nodes?
	ASSERT_EQ(grid->getNodes()->size(), (size_t)626);
	
	// Correct number of elements?
	ASSERT_EQ(grid->getElements()->size(), (size_t)1082);
	
	// Configure stream
	std::stringstream ss;
	ss.setf(ios::fixed);
	ss.precision(6);
	
	FileIO::OGSMeshIO::write (grid->getCFEMesh(), ss);
	
	compareToReference(ss.str().c_str(), QString("testMeshFromRaster_result.msh"));
}
