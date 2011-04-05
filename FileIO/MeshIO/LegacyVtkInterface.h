/**
 * \file LegacyVtkInterface.h
 * 05/04/2011 LB Initial implementation
 */

#ifndef LEGACYVTKINTERFACE_H
#define LEGACYVTKINTERFACE_H

#include <fstream>
#include <string>

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
	LegacyVtkInterface(Mesh_Group::CFEMesh* mesh, COutput* output);
	virtual ~LegacyVtkInterface();
	
	void WriteDataVTK(int number);

protected:
	void WriteVTKHeader(std::fstream&, int);
	void WriteVTKPointData(std::fstream&);
	void WriteVTKCellData(std::fstream&);
	void WriteVTKDataArrays(std::fstream&);
	void WriteELEVelocity(std::fstream &vtk_file);
	
	Mesh_Group::CFEMesh* _mesh;
	COutput* _output;
};

#endif // LEGACYVTKINTERFACE_H
