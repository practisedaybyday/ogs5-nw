/**
 * \file QtTestrunner.cpp
 * 2012-01-16 LB Initial implementation
 */

#include <QApplication>
#include <QtTest/QTest>
#include "Vtk/TestMeshFromRaster.h"

/// @brief If a function name is passed via the command line the reference is updated.
int main(int argc, char *argv[])
{	
	// Creating a non-gui (console) Qt application
	QApplication app(argc, argv, false);
	
	bool success = true;
	
	// Add your test here:
	TestMeshFromRaster testMeshFromRaster(argc);
	success = success && !QTest::qExec(&testMeshFromRaster, argc, argv);
	
	return !success;
}