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

class LegacyVtkInterface
{	
public:
	LegacyVtkInterface(Mesh_Group::CFEMesh* mesh, COutput* output);
	virtual ~LegacyVtkInterface();
	
	void WriteDataVTK(int number);

protected:
	void WriteVTKHeader(std::fstream&, int);
	void WriteVTKNodeData(std::fstream&);
	void WriteVTKElementData(std::fstream&);
	void WriteVTKValues(std::fstream&);
	void WriteELEVelocity(std::fstream &vtk_file);
	
	Mesh_Group::CFEMesh* _mesh;
	COutput* _output;
};

#endif // LEGACYVTKINTERFACE_H
