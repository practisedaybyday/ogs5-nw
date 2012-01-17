/**
 * \file QtTestrunner.cpp
 * 2012-01-16 LB Initial implementation
 */

#include <QApplication>
#include <QtTest/QTest>
#include "Vtk/TestMeshFromRaster.h"

int main(int argc, char *argv[])
{
	QApplication app(argc, argv);
	
	// Add your test here:
	TestMeshFromRaster testMeshFromRaster;
	QTest::qExec(&testMeshFromRaster, argc, argv);
	
	return 0;
}