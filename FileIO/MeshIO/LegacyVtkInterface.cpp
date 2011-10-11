/**
 * \file LegacyVtkInterface.cpp
 * 05/04/2011 LB Initial implementation
 *
 * Implementation of LegacyVtkInterface class
 */

// ** INCLUDES **
#include "LegacyVtkInterface.h"

#include "FEMEnums.h"
#include "ProcessInfo.h"
#include "fem_ele_std.h"
#include "matrix_class.h"
#include "msh_lib.h"
#include "msh_mesh.h"
#include "rf_mmp_new.h" // this is for class CMediumProperties, what else???
#include "rf_pcs.h"
#include "rf_pcs.h"
#ifdef GEM_REACT
#include "rf_REACT_GEM.h"
#endif // GEM_REACT

#include <string>
using namespace std;

LegacyVtkInterface::LegacyVtkInterface(MeshLib::CFEMesh* mesh,
                                       std::vector<std::string> pointArrayNames,
                                       std::vector<std::string> cellArrayNames,
                                       std::vector<std::string> materialPropertyArrayNames,
                                       std::string meshTypeName,
                                       ProcessInfo* processInfo)
	: _mesh(mesh),
	  _pointArrayNames(pointArrayNames),
	  _cellArrayNames(cellArrayNames),
	  _materialPropertyArrayNames(materialPropertyArrayNames),
	  _meshTypeName(meshTypeName),
	  _processInfo(processInfo)
{
	_processType = convertProcessTypeToString(processInfo->getProcessType());
	_mesh = FEMGet(_processType);
}

LegacyVtkInterface::~LegacyVtkInterface() {}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   04/2006 KG44 Implementation
   09/2006 KG44 Output for MPI - correct OUTPUT not yet implemented
   12/2008 NW Remove ios::app, Add PCS name to VTK file name
**************************************************************************/
void LegacyVtkInterface::WriteDataVTK(int number, double simulation_time,
                                      std::string baseFilename) const
{
#if defined(USE_MPI) || defined(USE_MPI_PARPROC) || defined(USE_MPI_REGSOIL)
	char tf_name[10];
	cout << "Process " << myrank << " in WriteDataVTK" << "\n";
#endif

	if(!_mesh)
	{
		cout << "Warning in LegacyVtkInterface::WriteVTKNodes - no MSH data" << endl;
		return;
	}

	if (_processType.compare("INVALID_PROCESS") != 0)
		baseFilename += "_" + _processType;

	stringstream ss;
	ss << number;
	baseFilename += ss.str();
	baseFilename += ".vtk";

	// LB if(!_new_file_opened) remove(baseFilename.c_str());
	fstream vtk_file (baseFilename.c_str(), ios::out);
	vtk_file.setf(ios::scientific,ios::floatfield);
	vtk_file.precision(12);
	if (!vtk_file.good())
	{
		cout << "Could not open file for writing: " << baseFilename << endl;
		return;
	}
	vtk_file.seekg(0L,ios::beg);
#ifdef SUPERCOMPUTER
// kg44 buffer the output
	char mybuffer [MY_IO_BUFSIZE * MY_IO_BUFSIZE];
	vtk_file.rdbuf()->pubsetbuf(mybuffer,MY_IO_BUFSIZE * MY_IO_BUFSIZE);
//
#endif

	this->WriteVTKHeader(vtk_file,number,simulation_time);
	this->WriteVTKPointData(vtk_file);
	this->WriteVTKCellData(vtk_file);
	this->WriteVTKDataArrays(vtk_file);

	vtk_file.close();
}

void LegacyVtkInterface::WriteVTKHeader(fstream &vtk_file,
                                        int time_step_number,
                                        double simulation_time) const
{
	vtk_file << "# vtk DataFile Version 3.0" << endl;
	vtk_file << "Unstructured Grid from OpenGeoSys" << endl;
	vtk_file << "ASCII"  << endl;
	vtk_file << "DATASET UNSTRUCTURED_GRID"  << endl;

	// time information
	// see http://www.visitusers.org/index.php?title=Time_and_Cycle_in_VTK_files
	vtk_file << "FIELD TimesAndCycles 2" << endl;
	vtk_file << "TIME 1 1 double" << endl;
	vtk_file << simulation_time << endl;
	vtk_file << "CYLCE 1 1 long" << endl;
	vtk_file << time_step_number << endl;
}

void LegacyVtkInterface::WriteVTKPointData(fstream &vtk_file) const
{
	const std::vector<MeshLib::CNode*> pointVector = _mesh->getNodeVector();
	vtk_file << "POINTS " << pointVector.size() << " double" << endl;

	const size_t n_pnts (pointVector.size());
	for(size_t i = 0; i < n_pnts; i++)
	{
		double const* const coords (pointVector[i]->getData());
		vtk_file << coords[0] << " " << coords[1] << " " << coords[2] << endl;
	}
}

void LegacyVtkInterface::WriteVTKCellData(fstream &vtk_file) const
{
	size_t numCells = _mesh->ele_vector.size();

	// count overall length of element vector
	long numAllPoints = 0;
	for(size_t i = 0; i < numCells; i++)
	{
		MeshLib::CElem* ele = _mesh->ele_vector[i];
		numAllPoints = numAllPoints + (ele->GetNodesNumber(false)) + 1;
	}

	// write elements
	vtk_file << "CELLS " << numCells << " " << numAllPoints << endl;
	for(size_t i = 0; i < numCells; i++)
	{
		MeshLib::CElem* ele = _mesh->ele_vector[i];

		// Write number of points per cell
		switch(ele->GetElementType())
		{
		case MshElemType::LINE:
			vtk_file << "2";
			break;
		case MshElemType::QUAD:
			vtk_file << "4";
			break;
		case MshElemType::HEXAHEDRON:
			vtk_file << "8";
			break;
		case MshElemType::TRIANGLE:
			vtk_file << "3";
			break;
		case MshElemType::TETRAHEDRON:
			vtk_file << "4";
			break;
		case MshElemType::PRISM:
			vtk_file << "6";
			break;
		case MshElemType::PYRAMID:
			vtk_file << "5";
			break;
		default:
			cerr << "COutput::WriteVTKElementData MshElemType not handled" << endl;
			break;
		}

		for(size_t j = 0; j < ele->GetNodesNumber(false); j++)
			vtk_file << " " << ele->nodes_index[j];

		vtk_file << endl;
	}
	vtk_file << endl;

	// write cell types
	vtk_file << "CELL_TYPES " << numCells << endl;

	for(size_t i = 0; i < numCells; i++)
	{
		MeshLib::CElem* ele = _mesh->ele_vector[i];

		// Write vtk cell type number (see vtkCellType.h)
		switch(ele->GetElementType())
		{
		case MshElemType::LINE:
			vtk_file << "3" << endl;
			break;
		case MshElemType::QUAD:
			vtk_file << "9" << endl;
			break;
		case MshElemType::HEXAHEDRON:
			vtk_file << "12" << endl;
			break;
		case MshElemType::TRIANGLE:
			vtk_file << "5" << endl;
			break;
		case MshElemType::TETRAHEDRON:
			vtk_file << "10" << endl;
			break;
		case MshElemType::PRISM: // VTK_WEDGE
			vtk_file << "13" << endl;
			break;
		case MshElemType::PYRAMID:
			vtk_file << "14" << endl;
			break;
		default:
			cerr << "COutput::WriteVTKElementData MshElemType not handled" << endl;
			break;
		}
	}
	vtk_file << endl;
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   04/2006 kg44 Implementation
   10/2006 WW Output secondary variables
   08/2008 OK MAT values
   06/2009 WW/OK WriteELEVelocity for different coordinate systems
**************************************************************************/
void LegacyVtkInterface::WriteVTKDataArrays(fstream &vtk_file) const
{
	std::vector<int> nod_value_index_vector(_pointArrayNames.size());
	long numNodes = _mesh->GetNodesNumber(false);

	// NODAL DATA
	vtk_file << "POINT_DATA " << numNodes << endl;
	for (size_t k = 0; k < _pointArrayNames.size(); k++)
	{
		string arrayName = _pointArrayNames[k];
		// Write x, y and z arrays as vectors
		if (arrayName.find("X") != string::npos)
		{
			vtk_file << "VECTORS " <<
			arrayName.substr(0, arrayName.size() - 2) << " double" << endl;
			string arrayNames[3];
			arrayNames[0] = arrayName;
			arrayNames[1] = arrayName.substr(0, arrayName.size() - 1).append("Y");
			arrayNames[2] = arrayName.substr(0, arrayName.size() - 1).append("Z");

			double vector3[3];
			for (long j = 0l; j < numNodes; j++)
			{
				for(int component = 0; component < 3; ++component)
				{
					CRFProcess* pcs = PCSGet(arrayNames[component], true);
					if (!pcs)
						continue;
					nod_value_index_vector[k] = pcs->GetNodeValueIndex(
					        arrayName);
					for (size_t i = 0; i < pcs->GetPrimaryVNumber(); i++)
						if (arrayName.compare(pcs->
						                      pcs_primary_function_name[i])
						    == 0)
						{
							nod_value_index_vector[k]++;
							break;
						}
					vector3[component] = pcs->GetNodeValue(
					        _mesh->nod_vector[j]->GetIndex(),
					        nod_value_index_vector
					        [k]);
				}
				vtk_file << vector3[0] << " " << vector3[1] << " " << vector3[2] <<
				endl;
			}
			k += 2;
		}
		else
		{
			CRFProcess* pcs = PCSGet(arrayName, true);
			if (!pcs)
				continue;
			nod_value_index_vector[k] = pcs->GetNodeValueIndex(arrayName);
			for (size_t i = 0; i < pcs->GetPrimaryVNumber(); i++)
				if (arrayName.compare(pcs->pcs_primary_function_name[i]) == 0)
				{
					nod_value_index_vector[k]++;
					break;
				}
			vtk_file << "SCALARS " << arrayName << " double 1" << endl;
			vtk_file << "LOOKUP_TABLE default" << endl;

			for (long j = 0l; j < numNodes; j++)
				if (nod_value_index_vector[k] > -1)
					vtk_file << pcs->GetNodeValue(
					        _mesh->nod_vector[j]->GetIndex(),
					        nod_value_index_vector[k])
					<< endl;
		}
	}
	//======================================================================
	// Saturation 2 for 1212 pp - scheme. 01.04.2009. WW
	// ---------------------------------------------------------------------
	CRFProcess* pcs = NULL;
	if (!_pointArrayNames.empty())                    //SB added
		pcs = PCSGet(_pointArrayNames[0], true);
	if (pcs && pcs->type == 1212)
	{
		size_t i = pcs->GetNodeValueIndex("SATURATION1");
		vtk_file << "SCALARS SATURATION2 double 1" << endl;
		vtk_file << "LOOKUP_TABLE default" << endl;
		for (long j = 0l; j < numNodes; j++)
		{
			double val_n = pcs->GetNodeValue(_mesh->nod_vector[j]->GetIndex(), i);
			vtk_file << 1.0 - val_n << endl;
		}
	}
//kg44 GEM node data
#ifdef GEM_REACT
	m_vec_GEM->WriteVTKGEMValues(vtk_file);    //kg44 export GEM internal variables like speciateion vector , phases etc
#endif
	// ELEMENT DATA
	// ---------------------------------------------------------------------
	bool wroteAnyEleData = false;              //NW
	if (!_cellArrayNames.empty())
	{
		CRFProcess* pcs = this->GetPCS_ELE(_cellArrayNames[0]);

		std::vector<int> ele_value_index_vector(_cellArrayNames.size());
		if (_cellArrayNames[0].size() > 0)
			for (size_t i = 0; i < _cellArrayNames.size(); i++)
				ele_value_index_vector[i] = pcs->GetElementValueIndex(
				        _cellArrayNames[i]);

		vtk_file << "CELL_DATA " << (long) _mesh->ele_vector.size() << endl;
		wroteAnyEleData = true;
		//....................................................................
		//
		for (size_t k = 0; k < _cellArrayNames.size(); k++)
		{
			// JTARON 2010, "VELOCITY" should only write as vector, scalars handled elswhere
			if (_cellArrayNames[k].compare("VELOCITY") == 0)
			{
				vtk_file << "VECTORS velocity double " << endl;
				this->WriteELEVelocity(vtk_file); //WW/OK
			}
			// PRINT CHANGING (OR CONSTANT) PERMEABILITY TENSOR?   // JTARON 2010
			else if (_cellArrayNames[k].compare("PERMEABILITY") == 0)
			{
				vtk_file << "VECTORS permeability double " << endl;
				for (int j = 0l; j < (long) _mesh->ele_vector.size(); j++)
				{
					MeshLib::CElem* ele = _mesh->ele_vector[j];
					CMediumProperties* MediaProp =
					        mmp_vector[ele->GetPatchIndex()];
					for (size_t i = 0; i < 3; i++)
						vtk_file <<
						MediaProp->PermeabilityTensor(j)[i * 3 + i] << " ";
					vtk_file << endl;
				}
			}
			else if (ele_value_index_vector[k] > -1)
			{
				// NOW REMAINING SCALAR DATA  // JTARON 2010, reconfig
				vtk_file << "SCALARS " << _cellArrayNames[k] << " double 1"
				         << endl;
				vtk_file << "LOOKUP_TABLE default" << endl;
				for (size_t i = 0; i < _mesh->ele_vector.size(); i++)
					vtk_file << pcs->GetElementValue(i,
					                                 ele_value_index_vector[k])
					         << endl;
			}
		}
		//--------------------------------------------------------------------
		ele_value_index_vector.clear();
	}
	//======================================================================
	// MAT data
	if (!_materialPropertyArrayNames.empty())
	{
		int mmp_id = -1;
		if (_materialPropertyArrayNames[0].compare("POROSITY") == 0)
			mmp_id = 0;
		// Let's say porosity
		// write header for cell data
		if (!wroteAnyEleData)
			vtk_file << "CELL_DATA " << _mesh->ele_vector.size() << endl;
		wroteAnyEleData = true;
		for (size_t i = 0; i < _mesh->ele_vector.size(); i++)
		{
			MeshLib::CElem* ele = _mesh->ele_vector[i];
			double mat_value = 0.0;
			switch (mmp_id)
			{
			case 0:
				mat_value = mmp_vector[ele->GetPatchIndex()]->Porosity(i, 0.0);
				break;
			default:
				cout << "COutput::WriteVTKValues: no MMP values specified" << endl;
				break;
			}
			vtk_file << mat_value << endl;
		}
	}
	// PCH: Material groups from .msh just for temparary purpose
	if (mmp_vector.size() > 1)
	{
		// write header for cell data
		if (!wroteAnyEleData)               //NW: check whether the header has been already written
			vtk_file << "CELL_DATA " << _mesh->ele_vector.size() << endl;
		wroteAnyEleData = true;

		vtk_file << "SCALARS " << "MatGroup" << " int 1" << endl;
		vtk_file << "LOOKUP_TABLE default" << endl;
		for (size_t i = 0; i < _mesh->ele_vector.size(); i++)
		{
			MeshLib::CElem* ele = _mesh->ele_vector[i];
			vtk_file << ele->GetPatchIndex() << endl;
		}
	}
}

/**************************************************************************
   FEMLib-Method:
   06/2009 WW/OK Implementation
**************************************************************************/
inline void LegacyVtkInterface::WriteELEVelocity(fstream &vtk_file) const
{
	int vel_ind[3] = {0, 1, 2};

	// 1D
	if(_mesh->GetCoordinateFlag() / 10 == 1)
	{
		// 0 y 0
		if(_mesh->GetCoordinateFlag() % 10 == 1)
		{
			vel_ind[0] = 1;
			vel_ind[1] = 0;
		}
		// 0 0 z
		else if(_mesh->GetCoordinateFlag() % 10 == 2)
		{
			vel_ind[0] = 2;
			vel_ind[2] = 0;
		}
	}
	// 2D
	if(_mesh->GetCoordinateFlag() / 10 == 2)
	{
		// 0 y z
		if(_mesh->GetCoordinateFlag() % 10 == 1)
		{
			vel_ind[0] = 1;
			vel_ind[1] = 2;
		}
		// x 0 z
		else if(_mesh->GetCoordinateFlag() % 10 == 2)
		{
			vel_ind[0] = 0;
			vel_ind[1] = 2;
			vel_ind[2] = 1;
		}
	}

	for(long i = 0; i < (long)_mesh->ele_vector.size(); i++)
	{
		for(int k = 0; k < 3; k++)
			vtk_file << ele_gp_value[i]->Velocity(vel_ind[k],0) << " ";
		vtk_file << endl;
	}
}

CRFProcess* LegacyVtkInterface::GetPCS_ELE(const string &var_name) const
{
	string pcs_var_name;
	CRFProcess* pcs = NULL;
	if (_processInfo->getProcessType() != INVALID_PROCESS)
		pcs = PCSGet(_processInfo->getProcessType());
	else if (_meshTypeName.size() > 0)
		pcs = PCSGet(_meshTypeName);

	if (!pcs)
		for (size_t i = 0; i < pcs_vector.size(); i++)
		{
			pcs = pcs_vector[i];
			for (int j = 0; j < pcs->pcs_number_of_evals; j++)
			{
				pcs_var_name = pcs->pcs_eval_name[j];
				if (pcs_var_name.compare(var_name) == 0)
					return pcs;
			}
		}
	return pcs;
}
