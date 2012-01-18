/**
 * \file QtTestrunner.cpp
 * 2012-01-16 LB Initial implementation
 */

#include <QApplication>
#include <QtTest/QTest>
#include "Vtk/TestMeshFromRaster.h"

int main(int argc, char *argv[])
{
	// Creating a non-gui (console) Qt application
	QApplication app(argc, argv, false);
	
	// Add your test here:
	TestMeshFromRaster testMeshFromRaster;
	QTest::qExec(&testMeshFromRaster, argc, argv);
	
	return 0;
}