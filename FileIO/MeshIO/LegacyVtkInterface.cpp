/**
 * \file LegacyVtkInterface.cpp
 * 05/04/2011 LB Initial implementation
 * 
 * Implementation of LegacyVtkInterface class
 */

// ** INCLUDES **
#include "LegacyVtkInterface.h"

#include "msh_mesh.h"
#include "msh_lib.h"
#include "Output.h"
#include "rf_pcs.h"
#include "rf_mmp_new.h" // this is for class CMediumProperties, what else???
#include "fem_ele_std.h"
#include "matrix_class.h"

#include <string>
using namespace std;

LegacyVtkInterface::LegacyVtkInterface(Mesh_Group::CFEMesh* mesh, COutput* output)
: _mesh(mesh), _output(output)
{
	
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
void LegacyVtkInterface::WriteDataVTK(int number)
{
	string filename = _output->file_base_name;
	char number_char[10];
	sprintf(number_char,"%i",number);
	string number_string = number_char;

#if defined(USE_MPI) || defined(USE_MPI_PARPROC) || defined(USE_MPI_REGSOIL)
	char tf_name[10];
	cout << "Process " << myrank << " in WriteDataVTK" << "\n";
#endif

	_mesh = FEMGet(convertProcessTypeToString (_output->getProcessType()));
	if(!_mesh)
	{
		cout << "Warning in COutput::WriteVTKNodes - no MSH data" << endl;
		return;
	}
//--------------------------------------------------------------------
// File handling
//  if(pcs_type_name.size()>0) // PCS
//    vtk_file_name += "_" + pcs_type_name;
	if(_output->getProcessType() != INVALID_PROCESS)        // PCS
		filename += "_" + convertProcessTypeToString (_output->getProcessType());

	filename += number_string ;

#if defined(USE_MPI) || defined(USE_MPI_PARPROC) || defined(USE_MPI_REGSOIL)
// kg44 removed "_0" as this make it impossible for visit/paraview to identify the cycle number
// sprintf(tf_name, "%d", myrank);
// filename += "_" + string(tf_name);
// std::cout << "VTK filename: " << filename << endl;
#endif
	filename += ".vtk";
												//KG44
// LB if(!_new_file_opened) remove(filename.c_str());
												//NW remove ios::app
	fstream vtk_file (filename.data(),ios::out);
	vtk_file.setf(ios::scientific,ios::floatfield);
	vtk_file.precision(12);
	if (!vtk_file.good()) return;
	vtk_file.seekg(0L,ios::beg);
#ifdef SUPERCOMPUTER
// kg44 buffer the output
	char mybuffer [MY_IO_BUFSIZE*MY_IO_BUFSIZE];
	vtk_file.rdbuf()->pubsetbuf(mybuffer,MY_IO_BUFSIZE*MY_IO_BUFSIZE);
//
#endif
//--------------------------------------------------------------------
	this->WriteVTKHeader(vtk_file,number);
	this->WriteVTKNodeData(vtk_file);
	this->WriteVTKElementData(vtk_file);
	this->WriteVTKValues(vtk_file);
//======================================================================
// vtk
	vtk_file.close();                              // kg44 close file
}


/**************************************************************************
FEMLib-Method:
Programing:
04/2006 KG44 Implementation
12/2009 KG44 added time information to header
**************************************************************************/
void LegacyVtkInterface::WriteVTKHeader(fstream &vtk_file, int time_step_number)
{
	
	vtk_file << "# vtk DataFile Version 3.0" << endl;
	vtk_file << "Unstructured Grid from OpenGeoSys " << OGS_VERSION << endl;
	vtk_file << "ASCII"  << endl;
	vtk_file << "DATASET UNSTRUCTURED_GRID"  << endl;
	
	// time information
	(void)time_step_number;
	//vtk_file << "FIELD TimesAndCycles 2" << endl;
	//vtk_file << "TIME 1 1 double" << endl;
	//vtk_file << _time << endl;
	//vtk_file << "CYLCE 1 1 long" << endl;
	//vtk_file << time_step_number << endl;

}


void LegacyVtkInterface::WriteVTKNodeData(fstream &vtk_file)
{
// header for node data
//   CFEMesh* _mesh = GetMSH(); //WW
//   _mesh = GetMSH(); //WW
	//KG44 12/2009 should be double !!!!
	vtk_file << "POINTS "<< _mesh->GetNodesNumber(false) << " double" << endl;

	for(long i=0; i < _mesh->GetNodesNumber(false) ;i++)
	{
		CNode* m_nod = _mesh->nod_vector[i];
		vtk_file << m_nod->X() << " " << m_nod->Y() << " " << m_nod->Z() << endl;
	}
}


void LegacyVtkInterface::WriteVTKElementData(fstream &vtk_file)
{

	int j;
	long no_all_elements =0;
	CElem* m_ele = NULL;
//  CFEMesh* _mesh = GetMSH(); //WW
//  _mesh = GetMSH(); //WW

	size_t numCells = _mesh->ele_vector.size();

// count overall length of element vector
	for(size_t i=0; i < numCells; i++)
	{
		m_ele = _mesh->ele_vector[i];
		no_all_elements=no_all_elements+(m_ele->GetNodesNumber(false))+1;
	}

// write element header
	vtk_file << "CELLS " << numCells << " " << no_all_elements << endl;

// write elements
	for(size_t i=0; i < numCells; i++)
	{
		m_ele = _mesh->ele_vector[i];

		switch(m_ele->GetElementType())
		{
			case MshElemType::LINE:                  // vtk_line=3
			vtk_file << "2" ;
			break;
			case MshElemType::QUAD:                  // quadrilateral=9
			vtk_file << "4";
			break;
			case MshElemType::HEXAHEDRON:            // hexahedron=12
			vtk_file << "8";
			break;
			case MshElemType::TRIANGLE:              // triangle=5
			vtk_file << "3";
			break;
			case MshElemType::TETRAHEDRON:           // tetrahedron=10
			vtk_file << "4";
			break;
			case MshElemType::PRISM:                 // wedge=13
			vtk_file << "6";
			break;
			default:
			std::cerr << "COutput::WriteVTKElementData MshElemType not handled" << std::endl;
		}

		for(j=0;j<m_ele->GetNodesNumber(false);j++)
		{
			vtk_file << " " << m_ele->nodes_index[j];
		}
		vtk_file << endl;
	}

	vtk_file << endl;

// write cell types

// write cell_types header
	vtk_file << "CELL_TYPES " << numCells << endl;

	for(size_t i=0; i < numCells; i++)
	{
		m_ele = _mesh->ele_vector[i];

		switch(m_ele->GetElementType())
		{
			case MshElemType::LINE:                  // vtk_line=3
			vtk_file << "3" << endl ;
			break;
			case MshElemType::QUAD:                  // quadrilateral=9
			vtk_file << "9" << endl;
			break;
			case MshElemType::HEXAHEDRON:            // hexahedron=12
			vtk_file << "12" << endl;
			break;
			case MshElemType::TRIANGLE:              // triangle=5
			vtk_file << "5" << endl;
			break;
			case MshElemType::TETRAHEDRON:           // tetrahedron=10
			vtk_file << "10" << endl;
			break;
			case MshElemType::PRISM:                 // wedge=13
			vtk_file << "13" << endl;
			break;
			default:
			std::cerr << "COutput::WriteVTKElementData MshElemType not handled" << std::endl;
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
void LegacyVtkInterface::WriteVTKValues(fstream &vtk_file)
{
	CRFProcess* pcs = NULL;
	const size_t num_nod_values(_output->_nod_value_vector.size());
	const size_t num_ele_values(_output->_ele_value_vector.size());
	std::vector<int> nod_value_index_vector(num_nod_values);
	std::vector<int> ele_value_index_vector(num_ele_values);
	double val_n = 0.;
	double *tensor = NULL;                         // JTARON 2010, added for permeability output
	long numNodes = _mesh->GetNodesNumber(false);

// NODAL DATA
	vtk_file << "POINT_DATA " << numNodes << endl;
//WW
	for (size_t k = 0; k < num_nod_values; k++)
	{
		string arrayName = _output->_nod_value_vector[k];
		if (arrayName.find("X"))
		{
			vtk_file << "VECTORS " << arrayName.substr(0, arrayName.size() - 2) << " double" << endl;
			string arrayNames[3];
			arrayNames[0] = arrayName;
			arrayNames[1] = arrayName.substr(0, arrayName.size() - 1).append("Y");
			arrayNames[2] = arrayName.substr(0, arrayName.size() - 1).append("Z");

			double vector3[3];
			for (long j = 0l; j < numNodes; j++)
			{
				for(int component = 0; component < 3; ++component)
				{
					pcs = PCSGet(arrayNames[component], true);
					if (!pcs)
						continue;
					nod_value_index_vector[k] = pcs->GetNodeValueIndex(arrayName);
				//   if(_nod_value_vector[k].find("SATURATION")!=string::npos)
				//   NodeIndex[k]++;
					for (size_t i = 0; i < pcs->GetPrimaryVNumber(); i++)
					{
						if (arrayName.compare(pcs->pcs_primary_function_name[i]) == 0)
						{
							nod_value_index_vector[k]++;
							break;
						}
					}
					vector3[component] = pcs->GetNodeValue(_mesh->nod_vector[j]->GetIndex(),
						nod_value_index_vector[k]);
				}
				vtk_file << vector3[0] << " " << vector3[1] << " " << vector3[2] << endl;
			}
			k += 2;
		}
		else
		{
			pcs = PCSGet(arrayName, true);
			if (!pcs)
				continue;
			nod_value_index_vector[k] = pcs->GetNodeValueIndex(arrayName);
		//   if(_nod_value_vector[k].find("SATURATION")!=string::npos)
		//   NodeIndex[k]++;
			for (size_t i = 0; i < pcs->GetPrimaryVNumber(); i++)
			{
				if (arrayName.compare(pcs->pcs_primary_function_name[i]) == 0)
				{
					nod_value_index_vector[k]++;
					break;
				}
			}
			vtk_file << "SCALARS " << arrayName << " double 1" << endl;
			vtk_file << "LOOKUP_TABLE default" << endl;
		//....................................................................
			for (long j = 0l; j < numNodes; j++)
			{
				if (nod_value_index_vector[k] > -1)
					vtk_file << pcs->GetNodeValue(_mesh->nod_vector[j]->GetIndex(),
					nod_value_index_vector[k])
					<< endl;
			}
		}
	}
//======================================================================
// Saturation 2 for 1212 pp - scheme. 01.04.2009. WW
// ---------------------------------------------------------------------
	if (num_nod_values > 0)                        //SB added
		pcs = PCSGet(_output->_nod_value_vector[0], true);
	if (pcs && pcs->type == 1212)
	{
		size_t i = pcs->GetNodeValueIndex("SATURATION1");
		vtk_file << "SCALARS SATURATION2 double 1" << endl;
	//
		vtk_file << "LOOKUP_TABLE default" << endl;
	//....................................................................
		for (long j = 0l; j < numNodes; j++)
		{
												//WW
			val_n = pcs->GetNodeValue(_mesh->nod_vector[j]->GetIndex(), i);
			vtk_file << 1. - val_n << endl;
		}
	}
//kg44 GEM node data
#ifdef GEM_REACT
	m_vec_GEM->WriteVTKGEMValues(vtk_file);        //kg44 export GEM internal variables like speciateion vector , phases etc
#endif
// ELEMENT DATA
// ---------------------------------------------------------------------
	bool wroteAnyEleData = false;                  //NW
	if (num_ele_values > 0)
	{
		pcs = _output->GetPCS_ELE(_output->_ele_value_vector[0]);
		_output->GetELEValuesIndexVector(ele_value_index_vector);
		vtk_file << "CELL_DATA " << (long) _mesh->ele_vector.size() << endl;
		wroteAnyEleData = true;
	//....................................................................
	//
		for (size_t k = 0; k < num_ele_values; k++)
		{
		//    JTARON 2010, "VELOCITY" should only write as vector, scalars handled elswhere
			if (_output->_ele_value_vector[k].compare("VELOCITY") == 0)
			{
				vtk_file << "VECTORS velocity double " << endl;
				this->WriteELEVelocity(vtk_file);           //WW/OK
			}
		//	  PRINT CHANGING (OR CONSTANT) PERMEABILITY TENSOR?   // JTARON 2010
			else if (_output->_ele_value_vector[k].compare("PERMEABILITY") == 0)
			{
				vtk_file << "VECTORS permeability double " << endl;
				CMediumProperties* MediaProp = NULL;
				CElem* m_ele = NULL;
				for (int j = 0l; j < (long) _mesh->ele_vector.size(); j++)
				{
					m_ele = _mesh->ele_vector[j];
					MediaProp = mmp_vector[m_ele->GetPatchIndex()];
					tensor = MediaProp->PermeabilityTensor(j);
					for (size_t i = 0; i < 3; i++)
						vtk_file << tensor[i * 3 + i] << " ";
					vtk_file << endl;
				}
			}
			else if (ele_value_index_vector[k] > -1)
			{
			//	  NOW REMAINING SCALAR DATA  // JTARON 2010, reconfig
				vtk_file << "SCALARS " << _output->_ele_value_vector[k] << " double 1"
					<< endl;
				vtk_file << "LOOKUP_TABLE default" << endl;
				for (size_t i = 0; i < _mesh->ele_vector.size(); i++)
					vtk_file << pcs->GetElementValue(i,
					ele_value_index_vector[k]) << endl;
			}
		}
	//--------------------------------------------------------------------
		ele_value_index_vector.clear();
	}
//======================================================================
// MAT data
	double mat_value = 0.0;                        //OK411
	CMediumProperties* m_mmp = NULL;
	CElem* m_ele = NULL;
	int mmp_id = -1;
	if (_output->mmp_value_vector.size() > 0)
	{
	// Identify MAT value
		if (_output->mmp_value_vector[0].compare("POROSITY") == 0)
			mmp_id = 0;
	// Let's say porosity
	// write header for cell data
		if (!wroteAnyEleData)
			vtk_file << "CELL_DATA " << _mesh->ele_vector.size() << endl;
		wroteAnyEleData = true;
		for (size_t i = 0; i < _mesh->ele_vector.size(); i++)
		{
			m_ele = _mesh->ele_vector[i];
			m_mmp = mmp_vector[m_ele->GetPatchIndex()];
			switch (mmp_id)
			{
				case 0:
				mat_value = m_mmp->Porosity(i, 0.0);
				break;
				default:
				cout << "COutput::WriteVTKValues: no MMP values specified"
					<< endl;
			}
			vtk_file << mat_value << endl;
		}
	}
// PCH: Material groups from .msh just for temparary purpose
	if (mmp_vector.size() > 1)
	{
	// write header for cell data
		if (!wroteAnyEleData)                       //NW: check whether the header has been already written
			vtk_file << "CELL_DATA " << _mesh->ele_vector.size() << endl;
		wroteAnyEleData = true;
	// header now scalar data
		vtk_file << "SCALARS " << "MatGroup" << " int 1" << endl;
		vtk_file << "LOOKUP_TABLE default" << endl;
		for (size_t i = 0; i < _mesh->ele_vector.size(); i++)
		{
			m_ele = _mesh->ele_vector[i];
			vtk_file << m_ele->GetPatchIndex() << endl;
		}
	}
}

/**************************************************************************
FEMLib-Method:
06/2009 WW/OK Implementation
**************************************************************************/
inline void LegacyVtkInterface::WriteELEVelocity(fstream &vtk_file)
{
	int k;
	int vel_ind[3];

	vel_ind[0] = 0;
	vel_ind[1] = 1;
	vel_ind[2] = 2;
// 1D
	if(_mesh->GetCoordinateFlag()/10==1)
	{
	// 0 y 0
		if(_mesh->GetCoordinateFlag()%10==1)
		{
			vel_ind[0] = 1;
			vel_ind[1] = 0;
		}
	// 0 0 z
		else if(_mesh->GetCoordinateFlag()%10==2)
		{
			vel_ind[0] = 2;
			vel_ind[2] = 0;
		}
	}
// 2D
	if(_mesh->GetCoordinateFlag()/10==2)
	{
	// 0 y z
		if(_mesh->GetCoordinateFlag()%10==1)
		{
			vel_ind[0] = 1;
			vel_ind[1] = 2;
		}
	// x 0 z
		else if(_mesh->GetCoordinateFlag()%10==2)
		{
			vel_ind[0] = 0;
			vel_ind[1] = 2;
			vel_ind[2] = 1;
		}
	}

	for(long i=0;i<(long)_mesh->ele_vector.size();i++)
	{
		for(k=0; k<3; k++)
			vtk_file << ele_gp_value[i]->Velocity(vel_ind[k],0) << " ";
		vtk_file << endl;
	}
}
