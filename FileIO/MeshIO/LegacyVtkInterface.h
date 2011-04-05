/**
 * \file LegacyVtkInterface.h
 * 05/04/2011 LB Initial implementation
 */

#ifndef LEGACYVTKINTERFACE_H
#define LEGACYVTKINTERFACE_H

#include <fstream>
#include <string>
#include <vector>

namespace Mesh_Group
{
   class CFEMesh;
}
class COutput;

/// @brief Writes a legacy ascii vtk file of a mesh.
// TODO decouple from COutput
class LegacyVtkInterface
{	
public:
	LegacyVtkInterface(Mesh_Group::CFEMesh* mesh, COutput* output,
		std::string processType,
		std::vector<std::string> pointArrayNames,
		std::vector<std::string> cellArrayNames,
		std::vector<std::string> materialPropertyArrayNames);
	virtual ~LegacyVtkInterface();
	
	void WriteDataVTK(int number, std::string baseFilename);

protected:
	void WriteVTKHeader(std::fstream&, int);
	void WriteVTKPointData(std::fstream&);
	void WriteVTKCellData(std::fstream&);
	void WriteVTKDataArrays(std::fstream&);
	void WriteELEVelocity(std::fstream &vtk_file);
	
	Mesh_Group::CFEMesh* _mesh;
	COutput* _output;
	std::string _processType;
	std::vector<std::string> _pointArrayNames;
	std::vector<std::string> _cellArrayNames;
	std::vector<std::string> _materialPropertyArrayNames;
};

#endif // LEGACYVTKINTERFACE_H
