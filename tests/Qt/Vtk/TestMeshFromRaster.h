/*
 * testMeshFromRaster.cpp
 * 2012/01/13 KR Initial implementation
 */

#ifndef TESTMESHFROMRASTER_H
#define TESTMESHFROMRASTER_H

#include "Configure.h"
#include "GridAdapter.h"
#include "VtkMeshConverter.h"
#include "VtkGeoImageSource.h"

#include <iostream>
#include <sstream>

#include <vtkSmartPointer.h>
#include <vtkImageData.h>

#include <QFile>
#include <QString>
#include <QTest>
#include <QTextStream>

#include "gdiff.h" // This should be the last include.

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

class TestMeshFromRaster : public QObject
{
	Q_OBJECT
	
private slots:
	void test()
	{
		QString fileName(SOURCEPATH);
		fileName += "/UTL/VTK/testMeshFromRaster.asc";
		vtkSmartPointer<VtkGeoImageSource> geo_image = vtkSmartPointer<VtkGeoImageSource>::New();
		geo_image->setImageFilename(fileName);
		vtkSmartPointer<vtkImageData> image = geo_image->GetOutput();

		GridAdapter* grid = VtkMeshConverter::convertImgToMesh(image, geo_image->getOrigin(),
			geo_image->getSpacing(), MshElemType::TRIANGLE, UseIntensityAs::ELEVATION);
		
		std::stringstream ss;
		FileIO::OGSMeshIO::write (grid->getCFEMesh(), ss);
		
		
		QString refFile(SOURCEPATH);
		refFile += "/UTL/VTK/testMeshFromRaster_result.msh";
		
		QFile qFile(refFile);
		if(!qFile.open(QIODevice::ReadOnly | QIODevice::Text))
			QFAIL(QString("Reference file %1 could not be read.").arg(refFile).toAscii());
		
		QString refFileContent = qFile.readAll();
		
		// File compare
		diff_match_patch gdiff;
		QString newContent = ss.str().c_str();
		QList<Diff> diffs = gdiff.diff_main(newContent, refFileContent);
		
		if(diffs.length() > 0)
		{
			QString htmlOutput = gdiff.diff_prettyHtml(diffs);
			QString htmlOutputFilename = "test.html";
			QFile htmlFile(htmlOutputFilename);
			if (!htmlFile.open(QIODevice::WriteOnly | QIODevice::Text))
				QFAIL("Html output file could not be written.");
			QTextStream htmlStream(&htmlFile);
			htmlStream << htmlOutput;
			QFAIL(QString("File compare failed. See %1 for the differences.")
				.arg(htmlOutputFilename).toAscii());
		}
	}
};

#endif // TESTMESHFROMRASTER_H