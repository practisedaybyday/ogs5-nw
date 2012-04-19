/**
 * \file FEM/Output.cpp
 * 05/04/2011 LB Refactoring: Moved from rf_out_new.h
 *
 * Implementation of Output class
 */

// ** INCLUDES **
#include "Configure.h"
#include "Output.h"

#include <fstream>
#include <iostream>
#include <string>

#include "Configure.h"

#include "FEMIO/GeoIO.h"
#include "GEOObjects.h"
#include "StringTools.h"
#include "fem_ele_std.h"
#include "files0.h"
#include "makros.h"
#include "mathlib.h"
#include "msh_lib.h"
#include "problem.h"
#include "rf_msp_new.h"
#include "rf_pcs.h"
#include "rf_pcs.h"
#include "rf_random_walk.h"
#include "rf_tim_new.h"
#include "vtk.h"

extern size_t max_dim;                            //OK411 todo

#ifdef CHEMAPP
#include "eqlink.h"
#endif

#if defined(USE_MPI) || defined(USE_MPI_PARPROC) || defined(USE_MPI_REGSOIL)
#include "par_ddc.h"
#endif
#ifdef SUPERCOMPUTER
// kg44 this is usefull for io-buffering as endl flushes the buffer
#define endl '\n'     // Introduced by WW. LB super bad programming style: this breaks platform independet IO
#define MY_IO_BUFSIZE 4096
#endif // SUPERCOMPUTER
#ifdef GEM_REACT
#include "rf_REACT_GEM.h"
#endif // GEM_REACT

using MeshLib::CFEMesh;
using namespace std;

COutput::COutput() :
	GeoInfo(GEOLIB::GEODOMAIN), ProcessInfo(), _id(0), out_amplifier(0.0),
	m_msh(NULL), nSteps(-1), _new_file_opened(false), dat_type_name("TECPLOT")
{
	tim_type_name = "TIMES";
	m_pcs = NULL;
	vtk = NULL; //NW
}

COutput::COutput(size_t id) :
	GeoInfo(GEOLIB::GEODOMAIN), ProcessInfo(), _id(id), out_amplifier(0.0),
	m_msh(NULL), nSteps(-1), _new_file_opened(false), dat_type_name("TECPLOT")
{
	tim_type_name = "TIMES";
	m_pcs = NULL;
	vtk = NULL; //NW
}

void COutput::init()
{
	if (getProcessType () == FiniteElement::INVALID_PROCESS)
	{
		std::cerr <<
		"COutput::init(): could not initialize process pointer (process type INVALID_PROCESS) and appropriate mesh"
		          << std::endl;
		std::cerr <<
		"COutput::init(): trying to fetch process pointer using msh_type_name ... " <<
		std::endl;
		std::cerr <<
		"COutput::init(): trying to fetch process pointer using msh_type_name ... " <<
		std::endl;
		if(msh_type_name.size() > 0)
		{
			_pcs = PCSGet(msh_type_name);
			if (_pcs)
				std::cerr << " successful" << std::endl;
			else
			{
				std::cerr << " failed" << std::endl;
				exit (1);
			}
		}
		else
			std::cerr << " failed" << std::endl;
	}

	m_msh = FEMGet(convertProcessTypeToString(getProcessType()));

    setInternalVarialbeNames(m_msh); //NW

}

COutput::~COutput()
{
	mmp_value_vector.clear();             //OK

	if (this->vtk != NULL)
		delete vtk;               //NW
}

const std::string& COutput::getGeoName () const
{
	return geo_name;
}

/**************************************************************************
   FEMLib-Method:
   Task: OUT read function
   Programing:
   06/2004 OK Implementation
   07/2004 WW Remove old files
   11/2004 OK string streaming by SB for lines
   03/2005 OK PCS_TYPE
   12/2005 OK DIS_TYPE
   12/2005 OK MSH_TYPE
   08/2008 OK MAT
   06/2010 TF formated, restructured, signature changed, use new GEOLIB data structures
   09/2010 TF signature changed, removed some variables
**************************************************************************/
ios::pos_type COutput::Read(std::ifstream& in_str,
                            const GEOLIB::GEOObjects& geo_obj, const std::string& unique_geo_name)
{
	std::string line_string;
	bool new_keyword = false;
	ios::pos_type position;
	bool new_subkeyword = false;
	std::string tec_file_name;
	ios::pos_type position_line;
	bool ok = true;
	std::stringstream in;
	string name;
	ios::pos_type position_subkeyword;

	// Schleife ueber alle Phasen bzw. Komponenten
	while (!new_keyword)
	{
		position = in_str.tellg();
		if (new_subkeyword)
			in_str.seekg(position_subkeyword, ios::beg);
		new_subkeyword = false;
		// SB new input		in_str.getline(buffer,MAX_ZEILE);
		// SB new input         line_string = buffer;
		line_string.clear();
		line_string = GetLineFromFile1(&in_str);
		if (line_string.size() < 1)
			break;

		if (Keyword(line_string))
			return position;

		// subkeyword found
		if (line_string.find("$NOD_VALUES") != string::npos)
		{
			while ((!new_keyword) && (!new_subkeyword))
			{
				position_subkeyword = in_str.tellg();
				//SB input with comments  in_str >> line_string>>ws;
				line_string = GetLineFromFile1(&in_str);
				if (line_string.find("#") != string::npos)
					return position;
				if (line_string.find("$") != string::npos)
				{
					new_subkeyword = true;
					break;
				}
				if (line_string.size() == 0)
					break;  //SB: empty line
				in.str(line_string);
				in >> name;
				//_alias_nod_value_vector.push_back(name);
                _nod_value_vector.push_back(name);
				in.clear();
			}

			continue;
		}
		//--------------------------------------------------------------------
		// subkeyword found //MX
		if (line_string.find("$PCON_VALUES") != string::npos)
		{
			while ((!new_keyword) && (!new_subkeyword))
			{
				position_subkeyword = in_str.tellg();
				in_str >> line_string >> ws;
				if (line_string.find("#") != string::npos)
					return position;
				if (line_string.find("$") != string::npos)
				{
					new_subkeyword = true;
					break;
				}
				if (line_string.size() == 0)
					break;
				_pcon_value_vector.push_back(line_string);
			}
			continue;
		}

		//--------------------------------------------------------------------
		// subkeyword found
		if (line_string.find("$ELE_VALUES") != string::npos)
		{
			ok = true;
			while (ok)
			{
				position_line = in_str.tellg();
				in_str >> line_string;
				if (SubKeyword(line_string))
				{
					in_str.seekg(position_line, ios::beg);
					ok = false;
					continue;
				}
				if (Keyword(line_string))
					return position;
				_ele_value_vector.push_back(line_string);
				in_str.ignore(MAX_ZEILE, '\n');
			}
			/*
			   // Commented by WW
			   // Remove files
			   tec_file_name = file_base_name + "_domain_ele" + TEC_FILE_EXTENSION;
			   remove(tec_file_name.c_str());
			 */
			continue;
		}
		//-------------------------------------------------------------------- // Added 03.2010 JTARON
		// subkeyword found
		if (line_string.find("$RWPT_VALUES") != string::npos)
		{
			while ((!new_keyword) && (!new_subkeyword))
			{
				position_subkeyword = in_str.tellg();
				line_string = GetLineFromFile1(&in_str);
				if (line_string.find("#") != string::npos)
					return position;
				if (line_string.find("$") != string::npos)
				{
					new_subkeyword = true;
					break;
				}
				if (line_string.size() == 0)
					break;  //SB: empty line
				in.str(line_string);
				in >> name;
				_rwpt_value_vector.push_back(name);
				in.clear();
			}
			continue;
		}

		//subkeyword found
		if (line_string.find("$GEO_TYPE") != string::npos)
		{
			FileIO::GeoIO::readGeoInfo (this,
			                            in_str,
			                            geo_name,
			                            geo_obj,
			                            unique_geo_name);
			continue;
		}

		// subkeyword found
		if (line_string.find("$TIM_TYPE") != string::npos)
		{
			while ((!new_keyword) && (!new_subkeyword))
			{
				position_subkeyword = in_str.tellg();
				in_str >> line_string;
				if (line_string.size() == 0) //SB
					break;
				if (line_string.find("#") != string::npos)
				{
					new_keyword = true;
					break;
				}
				if (line_string.find("$") != string::npos)
				{
					new_subkeyword = true;
					break;
				}
				if (line_string.find("STEPS") != string::npos)
				{
					in_str >> nSteps;
					tim_type_name = "STEPS"; //OK
					break; //kg44 I guess that was missing..otherwise it pushes back a time_vector!
				}
				// JT 2010, reconfigured (and added RWPT)... didn't work
				if (line_string.find("STEPPING") != string::npos)
				{
					double stepping_length, stepping_end, stepping_current;
					in_str >> stepping_length >> stepping_end;
					stepping_current = stepping_length;
					while (stepping_current <= stepping_end)
					{
						time_vector.push_back(stepping_current);
						//						rwpt_time_vector.push_back(stepping_current);
						stepping_current += stepping_length;
					}
				}
				else
					time_vector.push_back(strtod(line_string.data(), NULL));
				//					rwpt_time_vector.push_back(strtod(line_string.data(), NULL));
				in_str.ignore(MAX_ZEILE, '\n');
			}
			continue;
		}

		// subkeyword found
		if (line_string.find("$DAT_TYPE") != string::npos)
		{
			in_str >> dat_type_name;
			in_str.ignore(MAX_ZEILE, '\n');
			continue;
		}

		// subkeyword found
		if (line_string.find("$AMPLIFIER") != string::npos)
		{
			in_str >> out_amplifier;
			in_str.ignore(MAX_ZEILE, '\n');
			continue;
		}

		// subkeyword found
		if (line_string.find("$PCS_TYPE") != string::npos)
		{
			std::string tmp_pcs_type_name;
			in_str >> tmp_pcs_type_name;
			setProcessType(FiniteElement::convertProcessType(tmp_pcs_type_name));
			in_str.ignore(MAX_ZEILE, '\n');
			/* // Comment by WW
			   // Remove files
			   tec_file_name = pcs_type_name + "_" + "domain" + "_line" + TEC_FILE_EXTENSION;
			   remove(tec_file_name.c_str());
			   tec_file_name = pcs_type_name + "_" + "domain" + "_tri" + TEC_FILE_EXTENSION;
			   remove(tec_file_name.c_str());
			   tec_file_name = pcs_type_name + "_" + "domain" + "_quad" + TEC_FILE_EXTENSION;
			   remove(tec_file_name.c_str());
			   tec_file_name = pcs_type_name + "_" + "domain" + "_tet" + TEC_FILE_EXTENSION;
			   remove(tec_file_name.c_str());
			   tec_file_name = pcs_type_name + "_" + "domain" + "_pris" + TEC_FILE_EXTENSION;
			   remove(tec_file_name.c_str());
			   tec_file_name = pcs_type_name + "_" + "domain" + "_hex" + TEC_FILE_EXTENSION;
			   remove(tec_file_name.c_str());
			 */
			continue;
		}

		// subkeyword found
		if (line_string.find("$DIS_TYPE") != string::npos)
		{
			std::string dis_type_name;
			in_str >> dis_type_name;
			setProcessDistributionType (FiniteElement::convertDisType(dis_type_name));
			in_str.ignore(MAX_ZEILE, '\n');
			continue;
		}

		// subkeyword found
		if (line_string.find("$MSH_TYPE") != string::npos)
		{
			in_str >> msh_type_name;
			in_str.ignore(MAX_ZEILE, '\n');
			continue;
		}

		//OK
		if (line_string.find("$MMP_VALUES") != string::npos)
		{
			ok = true;
			while (ok)
			{
				position_line = in_str.tellg();
				in_str >> line_string;
				if (SubKeyword(line_string))
				{
					in_str.seekg(position_line, ios::beg);
					ok = false;
					continue;
				}
				if (Keyword(line_string))
					return position;
				mmp_value_vector.push_back(line_string);
				in_str.ignore(MAX_ZEILE, '\n');
			}
			continue;
		}

		//OK
		if (line_string.find("$MFP_VALUES") != string::npos)
		{
			ok = true;
			while (ok)
			{
				position_line = in_str.tellg();
				in_str >> line_string;
				if (SubKeyword(line_string))
				{
					in_str.seekg(position_line, ios::beg);
					ok = false;
					continue;
				}
				if (Keyword(line_string))
					return position;
				mfp_value_vector.push_back(line_string);
				in_str.ignore(MAX_ZEILE, '\n');
			}

			continue;
		}
	}
	return position;
}

/**************************************************************************
   FEMLib-Method:
   Task: write function
   Programing:
   06/2004 OK Implementation
   01/2005 OK Extensions
   12/2005 OK DIS_TYPE
   12/2005 OK MSH_TYPE
   12/2008 NW DAT_TYPE
   05/2009 OK bug fix STEPS
**************************************************************************/
void COutput::Write(fstream* out_file)
{
	//--------------------------------------------------------------------
	// KEYWORD
	*out_file << "#OUTPUT" << endl;
	//--------------------------------------------------------------------
	// PCS_TYPE
	*out_file << " $PCS_TYPE" << std::endl << "  ";
	*out_file << convertProcessTypeToString(getProcessType()) << std::endl;
	//--------------------------------------------------------------------
	// NOD_VALUES
	*out_file << " $NOD_VALUES" << endl;
	size_t nod_value_vector_size(_nod_value_vector.size());
	for (size_t i = 0; i < nod_value_vector_size; i++)
		*out_file << "  " << _nod_value_vector[i] << endl;
	//--------------------------------------------------------------------
	// ELE_VALUES
	*out_file << " $ELE_VALUES" << endl;
	size_t ele_value_vector_size (_ele_value_vector.size());
	for (size_t i = 0; i < ele_value_vector_size; i++)
		*out_file << "  " << _ele_value_vector[i] << endl;
	//--------------------------------------------------------------------
	// GEO_TYPE
	*out_file << " $GEO_TYPE" << endl;
	*out_file << "  ";
	*out_file << getGeoTypeAsString() << " " << geo_name << endl;
	//--------------------------------------------------------------------
	// TIM_TYPE
	*out_file << " $TIM_TYPE" << endl;
	if (tim_type_name == "STEPS")
		*out_file << "  " << tim_type_name << " " << nSteps << endl;
	else
	{
		size_t time_vector_size (time_vector.size());
		for (size_t i = 0; i < time_vector_size; i++)
			*out_file << "  " << time_vector[i] << endl;
	}

	// DIS_TYPE
	//	if (_dis_type_name.size() > 0) {
	//		*out_file << " $DIS_TYPE" << endl;
	//		*out_file << "  ";
	//		*out_file << _dis_type_name << endl;
	//	}
	if (getProcessDistributionType() != FiniteElement::INVALID_DIS_TYPE)
	{
		*out_file << " $DIS_TYPE" << endl;
		*out_file << "  ";
		*out_file << convertDisTypeToString (getProcessDistributionType()) << std::endl;
	}

	// MSH_TYPE
	if (msh_type_name.size() > 0)
	{
		*out_file << " $MSH_TYPE" << endl;
		*out_file << "  ";
		*out_file << msh_type_name << endl;
	}
	//--------------------------------------------------------------------
	// DAT_TYPE
	*out_file << " $DAT_TYPE" << endl;
	*out_file << "  ";
	*out_file << dat_type_name << endl;
	//--------------------------------------------------------------------
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   08/2004 OK Implementation
   03/2005 OK MultiMSH
   08/2005 WW Changes for MultiMSH
   12/2005 OK VAR,MSH,PCS concept
   07/2007 NW Multi Mesh Type
**************************************************************************/
void COutput::NODWriteDOMDataTEC()
{
	int te = 0;
	string eleType;
	string tec_file_name;
#if defined(USE_MPI) || defined(USE_MPI_PARPROC) || defined(USE_MPI_REGSOIL)
	char tf_name[10];
	std::cout << "Process " << myrank << " in WriteDOMDataTEC" << "\n";
#endif
	//----------------------------------------------------------------------
	// Tests
	//OK4704
	if((_nod_value_vector.size() == 0) && (mfp_value_vector.size() == 0))
		return;
	//......................................................................
	// MSH
	//m_msh = FEMGet(pcs_type_name);
	//  m_msh = GetMSH();
	if(!m_msh)
	{
		cout << "Warning in COutput::NODWriteDOMDataTEC() - no MSH data" << endl;
		return;
	}
	//======================================================================
	vector<int> mesh_type_list;           //NW
	if (m_msh->getNumberOfLines() > 0)
		mesh_type_list.push_back(1);
	if (m_msh->getNumberOfQuads() > 0)
		mesh_type_list.push_back(2);
	if (m_msh->getNumberOfHexs () > 0)
		mesh_type_list.push_back(3);
	if (m_msh->getNumberOfTris () > 0)
		mesh_type_list.push_back(4);
	if (m_msh->getNumberOfTets () > 0)
		mesh_type_list.push_back(5);
	if (m_msh->getNumberOfPrisms () > 0)
		mesh_type_list.push_back(6);
	if (m_msh->getNumberOfPyramids() > 0)
		mesh_type_list.push_back(7);

	// Output files for each mesh type
	//NW
	for (int i = 0; i < (int)mesh_type_list.size(); i++)
	{
		te = mesh_type_list[i];
		//----------------------------------------------------------------------
		// File name handling
		tec_file_name = file_base_name + "_" + "domain";
		if(msh_type_name.size() > 0) // MultiMSH
			tec_file_name += "_" + msh_type_name;
		if(getProcessType() != FiniteElement::INVALID_PROCESS) // PCS
			tec_file_name += "_" + convertProcessTypeToString(getProcessType());
		//======================================================================
		switch (te)               //NW
		{
		case 1:
			tec_file_name += "_line";
			eleType = "QUADRILATERAL";
			break;
		case 2:
			tec_file_name += "_quad";
			eleType = "QUADRILATERAL";
			break;
		case 3:
			tec_file_name += "_hex";
			eleType = "BRICK";
			break;
		case 4:
			tec_file_name += "_tri";
			eleType = "QUADRILATERAL";
			break;
		case 5:
			tec_file_name += "_tet";
			eleType = "TETRAHEDRON";
			break;
		case 6:
			tec_file_name += "_pris";
			eleType = "BRICK";
			break;
		case 7:
			tec_file_name += "_pyra";
			eleType = "BRICK";
			break;
		}
		/*
		   if(m_msh->msh_no_line>0)
		   {
		      tec_file_name += "_line";
		      eleType = "QUADRILATERAL";
		     te=1;
		   }
		   else if (m_msh->msh_no_quad>0)
		   {
		      tec_file_name += "_quad";
		      eleType = "QUADRILATERAL";
		   te=2;
		   }
		   else if (m_msh->msh_no_hexs>0)
		   {
		   tec_file_name += "_hex";
		   eleType = "BRICK";
		   te=3;
		   }
		   else if (m_msh->msh_no_tris>0)
		   {
		   tec_file_name += "_tri";
		   //???Who was this eleType = "TRIANGLE";
		   eleType = "QUADRILATERAL";
		   te=4;
		   }
		   else if (m_msh->msh_no_tets>0)
		   {
		   tec_file_name += "_tet";
		   eleType = "TETRAHEDRON";
		   te=5;
		   }
		   else if (m_msh->msh_no_pris>0)
		   {
		   tec_file_name += "_pris";
		   eleType = "BRICK";
		   te=6;
		   }
		 */
#if defined(USE_MPI) || defined(USE_MPI_PARPROC) || defined(USE_MPI_REGSOIL)
		sprintf(tf_name, "%d", myrank);
		tec_file_name += "_" + string(tf_name);
		std::cout << "Tecplot filename: " << tec_file_name << endl;
#endif
		tec_file_name += TEC_FILE_EXTENSION;
		//WW
		if(!_new_file_opened)
			remove(tec_file_name.c_str());
		fstream tec_file (tec_file_name.data(),ios::app | ios::out);
		tec_file.setf(ios::scientific,ios::floatfield);
		tec_file.precision(12);
		if (!tec_file.good())
			return;
#ifdef SUPERCOMPUTER
		// kg44 buffer the output
		char mybuf1 [MY_IO_BUFSIZE * MY_IO_BUFSIZE];
		tec_file.rdbuf()->pubsetbuf(mybuf1,MY_IO_BUFSIZE * MY_IO_BUFSIZE);
#endif
		//
		WriteTECHeader(tec_file,te,eleType);
		WriteTECNodeData(tec_file);
		WriteTECElementData(tec_file,te);
		tec_file.close();         // kg44 close file
		//--------------------------------------------------------------------
		// tri elements
		// ***** 07/2010 TF commented out block since the global variable is always zero
		//    if(msh_no_tris>0){
		//    //string tec_file_name = pcs_type_name + "_" + "domain" + "_tri" + TEC_FILE_EXTENSION;
		//#ifdef SUPERCOMPUTER
		//// buffer the output
		//      char sxbuf1[MY_IO_BUFSIZE*MY_IO_BUFSIZE];
		//#endif
		//      string tec_file_name = file_base_name + "_" + "domain" + "_tri" + TEC_FILE_EXTENSION;
		//      fstream tec_file1 (tec_file_name.data(),ios::app|ios::out);
		//      tec_file1.setf(ios::scientific,ios::floatfield);
		//      tec_file1.precision(12);
		//      if (!tec_file1.good()) return;
		//#ifdef SUPERCOMPUTER
		//      tec_file1.rdbuf()->pubsetbuf(sxbuf1,MY_IO_BUFSIZE*MY_IO_BUFSIZE);
		//#endif
		//      //OK  tec_file1.clear();
		//      //OK  tec_file1.seekg(0L,ios::beg);
		//      WriteTECHeader(tec_file1,4,"TRIANGLE");
		//      WriteTECNodeData(tec_file1);
		//      WriteTECElementData(tec_file1,4);
		//      tec_file1.close(); // kg44 close file
		//    }
		//--------------------------------------------------------------------
		// quad elements
		// ***** 07/2010 TF commented out block since the global variable is always zero
		//    if(msh_no_quad>0){
		//      //string tec_file_name = pcs_type_name + "_" + "domain" + "_quad" + TEC_FILE_EXTENSION;
		//#ifdef SUPERCOMPUTER
		//      char sxbuf2[MY_IO_BUFSIZE*MY_IO_BUFSIZE];
		//#endif
		//      string tec_file_name = file_base_name + "_" + "domain" + "_quad" + TEC_FILE_EXTENSION;
		//      fstream tec_file (tec_file_name.data(),ios::app|ios::out);
		//
		//      tec_file.setf(ios::scientific,ios::floatfield);
		//      tec_file.precision(12);
		//      if (!tec_file.good()) return;
		//#ifdef SUPERCOMPUTER
		//      tec_file.rdbuf()->pubsetbuf(sxbuf2,MY_IO_BUFSIZE*MY_IO_BUFSIZE);
		//#endif
		//      WriteTECHeader(tec_file,2,"QUADRILATERAL");
		//      WriteTECNodeData(tec_file);
		//      WriteTECElementData(tec_file,2);
		//      tec_file.close(); // kg44 close file
		//    }
		//--------------------------------------------------------------------
		// tet elements
		// ***** 07/2010 TF commented out block since the global variable is always zero
		//    if(msh_no_tets>0){
		//      //string tec_file_name = pcs_type_name + "_" + "domain" + "_tet" + TEC_FILE_EXTENSION;
		//#ifdef SUPERCOMPUTER
		//      char sxbuf3[MY_IO_BUFSIZE*MY_IO_BUFSIZE];
		//#endif
		//
		//      string tec_file_name = file_base_name + "_" + "domain" + "_tet";
		//
		//#if defined(USE_MPI) || defined(USE_MPI_PARPROC) || defined(USE_MPI_REGSOIL)
		//      sprintf(tf_name, "%d", myrank);
		//      tec_file_name += "_" + string(tf_name);
		//#endif
		//
		//      tec_file_name += TEC_FILE_EXTENSION;
		//
		//      fstream tec_file (tec_file_name.data(),ios::app|ios::out);
		//
		//      tec_file.setf(ios::scientific,ios::floatfield);
		//      tec_file.precision(12);
		//      if (!tec_file.good()) return;
		//#ifdef SUPERCOMPUTER
		//      tec_file.rdbuf()->pubsetbuf(sxbuf3,MY_IO_BUFSIZE*MY_IO_BUFSIZE);
		//#endif
		//
		//      WriteTECHeader(tec_file,5,"TETRAHEDRON");
		//      WriteTECNodeData(tec_file);
		//      WriteTECElementData(tec_file,5);
		//      tec_file.close(); // kg44 close file
		//    }
		//--------------------------------------------------------------------
		//    // pris elements
		// ***** 07/2010 TF commented out block since the global variable is always zero
		//    if(msh_no_pris>0){
		//      //string tec_file_name = pcs_type_name + "_" + "domain" + "_pris" + TEC_FILE_EXTENSION;
		//#ifdef SUPERCOMPUTER
		//        char sxbuf4[MY_IO_BUFSIZE*MY_IO_BUFSIZE];
		//#endif
		//      string tec_file_name = file_base_name + "_" + "domain" + "_pris" + TEC_FILE_EXTENSION;
		//      fstream tec_file (tec_file_name.data(),ios::app|ios::out);
		//
		//      tec_file.setf(ios::scientific,ios::floatfield);
		//      tec_file.precision(12);
		//      if (!tec_file.good()) return;
		//#ifdef SUPERCOMPUTER
		//      tec_file.rdbuf()->pubsetbuf(sxbuf4,MY_IO_BUFSIZE*MY_IO_BUFSIZE);
		//#endif
		//
		//      WriteTECHeader(tec_file,6,"BRICK");
		//      WriteTECNodeData(tec_file);
		//      WriteTECElementData(tec_file,6);
		//      tec_file.close(); // kg44 close file
		//    }
		//--------------------------------------------------------------------
		// hex elements
		// ***** 07/2010 TF commented out block since the global variable is always zero
		//    if(msh_no_hexs>0){
		//      //string tec_file_name = pcs_type_name + "_" + "domain" + "_hex" + TEC_FILE_EXTENSION;
		//#ifdef SUPERCOMPUTER
		//        char sxbuf5[MY_IO_BUFSIZE*MY_IO_BUFSIZE];
		//#endif
		//
		//      string tec_file_name = file_base_name + "_" + "domain" + "_hex" + TEC_FILE_EXTENSION;
		//      fstream tec_file (tec_file_name.data(),ios::app|ios::out);
		//
		//
		//      tec_file.setf(ios::scientific,ios::floatfield);
		//      tec_file.precision(12);
		//      if (!tec_file.good()) return;
		//#ifdef SUPERCOMPUTER
		//      tec_file.rdbuf()->pubsetbuf(sxbuf5,MY_IO_BUFSIZE*MY_IO_BUFSIZE);
		//#endif
		//      WriteTECHeader(tec_file,3,"BRICK");
		//      WriteTECNodeData(tec_file);
		//      WriteTECElementData(tec_file,3);
		//      tec_file.close(); // kg44 close file
		//    }
	}
}

/**************************************************************************
   FEMLib-Method:
   Programing:
   08/2004 OK Implementation
   08/2004 WW Output node variables by their names given in .out file
   03/2005 OK MultiMSH
   08/2005 WW Correction of node index
   12/2005 OK Mass transport specifics
   OK ??? too many specifics
**************************************************************************/
void COutput::WriteTECNodeData(fstream &tec_file)
{
	const size_t nName(_nod_value_vector.size());
	double val_n = 0.;                    //WW
	int nidx, nidx_dm[3];
	vector<int> NodeIndex(nName);
	string nod_value_name;                //OK
	int timelevel;
	//	m_msh = GetMSH();
	CRFProcess* m_pcs_out = NULL;
	// MSH
	for (size_t k = 0; k < nName; k++)
	{
		m_pcs = PCSGet(_nod_value_vector[k], true);
		if (m_pcs != NULL)
		{
			NodeIndex[k] = m_pcs->GetNodeValueIndex(_nod_value_vector[k],true); // JT Latest.
		}
	}

	if (M_Process || MH_Process)          //WW
	{
		m_pcs = PCSGet("DISPLACEMENT_X1", true);
		nidx_dm[0] = m_pcs->GetNodeValueIndex("DISPLACEMENT_X1") + 1;
		nidx_dm[1] = m_pcs->GetNodeValueIndex("DISPLACEMENT_Y1") + 1;
		if (max_dim > 1)
			nidx_dm[2] = m_pcs->GetNodeValueIndex("DISPLACEMENT_Z1") + 1;
		else
			nidx_dm[2] = -1;
	}
	for (size_t j = 0; j < m_msh->GetNodesNumber(false); j++)
	{
		// XYZ
		double x[3] =
		{(m_msh->nod_vector[j]->getData())[0], (m_msh->nod_vector[j]->getData())[1],
		 (m_msh->nod_vector[j]->getData())[2]};
//      x[0] = m_msh->nod_vector[j]->X();
//      x[1] = m_msh->nod_vector[j]->Y();
//      x[2] = m_msh->nod_vector[j]->Z();
		// Amplifying DISPLACEMENTs
		if (M_Process || MH_Process) //WW

			for (size_t k = 0; k < max_dim + 1; k++)
				x[k] += out_amplifier * m_pcs->GetNodeValue(
				        m_msh->nod_vector[j]->GetIndex(), nidx_dm[k]);
		for (size_t i = 0; i < 3; i++)
			tec_file << x[i] << " ";

		// NOD values
		// Mass transport
		//     if(pcs_type_name.compare("MASS_TRANSPORT")==0){
		if (getProcessType() == FiniteElement::MASS_TRANSPORT)
			for (size_t i = 0; i < _nod_value_vector.size(); i++)
			{
				std::string nod_value_name = _nod_value_vector[i];
				for (size_t l = 0; l < pcs_vector.size(); l++)
				{
					m_pcs = pcs_vector[l];
					//					if (m_pcs->pcs_type_name.compare("MASS_TRANSPORT") == 0) {
					if (m_pcs->getProcessType () == FiniteElement::MASS_TRANSPORT)
					{
						timelevel = 0;
						for (size_t m = 0;
						     m < m_pcs->nod_val_name_vector.size(); m++)
							if (m_pcs->nod_val_name_vector[m].compare(
							            nod_value_name) == 0)
							{
								m_pcs_out = PCSGet(FiniteElement::MASS_TRANSPORT,
								                   nod_value_name);
								if (!m_pcs_out)
									continue;
								if (timelevel == 1)
								{
									nidx =
									        m_pcs_out->
									        GetNodeValueIndex(
									                nod_value_name)
									        +
									        timelevel;
									tec_file <<
									m_pcs_out->GetNodeValue(
									        m_msh->nod_vector[j
									        ]->GetIndex(),
									        nidx) << " ";
								}
								timelevel++;
							}
					}
				}
			}
		else
		{
			for (size_t k = 0; k < nName; k++)
			{
				m_pcs = GetPCS(_nod_value_vector[k]);
				if (m_pcs != NULL) //WW

					if (NodeIndex[k] > -1)
					{
						val_n = m_pcs->GetNodeValue(
						        //WW
						        m_msh->nod_vector[j]->GetIndex(),
						        NodeIndex[k]);
						tec_file << val_n << " ";
						if ((m_pcs->type == 1212 ||
						     m_pcs->type == 42) &&
						    _nod_value_vector[k].find(
						            //WW
						            "SATURATION") != string::npos)
							//WW
							tec_file << 1. - val_n << " ";
					}
			}
			//OK4704
			for (size_t k = 0; k < mfp_value_vector.size(); k++)
				//tec_file << MFPGetNodeValue(m_msh->nod_vector[j]->GetIndex(),mfp_value_vector[k]) << " "; //NB
				tec_file << MFPGetNodeValue(m_msh->nod_vector[j]->GetIndex(),
				                            mfp_value_vector[k],
				                            atoi(&mfp_value_vector[k][
				                                         mfp_value_vector[k].size()
				                                         - 1]) - 1) << " ";  //NB: MFP output for all phases
		}
		tec_file << endl;
	}
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   08/2004 OK Implementation
   03/2005 OK MultiMSH
   08/2005 WW Wite for MultiMSH
   12/2005 OK GetMSH
   07/2007 NW Multi Mesh Type
**************************************************************************/
void COutput::WriteTECElementData(fstream &tec_file,int e_type)
{
	for (size_t i = 0; i < m_msh->ele_vector.size(); i++)
	{
		if (!m_msh->ele_vector[i]->GetMark())
			continue;
		//NW
		if (m_msh->ele_vector[i]->GetElementType() == e_type)
			m_msh->ele_vector[i]->WriteIndex_TEC(tec_file);
	}
}

/**************************************************************************
   FEMLib-Method:
   Programing:
   08/2004 OK Implementation
   08/2004 WW Header by the names gives in .out file
   03/2005 OK MultiMSH
   04/2005 WW Output active elements only
   08/2005 WW Output by MSH
   12/2005 OK GetMSH
**************************************************************************/
void COutput::WriteTECHeader(fstream &tec_file,int e_type, string e_type_name)
{
	// MSH
	//	m_msh = GetMSH();

	//OK411
	size_t no_elements = 0;
	const size_t mesh_ele_vector_size (m_msh->ele_vector.size());
	for (size_t i = 0; i < mesh_ele_vector_size; i++)
		if (m_msh->ele_vector[i]->GetMark())
			if (m_msh->ele_vector[i]->GetElementType() == e_type)
				no_elements++;
	//--------------------------------------------------------------------
	// Write Header I: variables
	CRFProcess* pcs = NULL;               //WW
	const size_t nName (_nod_value_vector.size());
	tec_file << "VARIABLES  = \"X\",\"Y\",\"Z\"";
	for (size_t k = 0; k < nName; k++)
	{
		tec_file << ", \"" << _nod_value_vector[k] << "\"";
		//-------------------------------------WW
		pcs = GetPCS(_nod_value_vector[k]);
		if (pcs != NULL)
			if ((pcs->type == 1212 ||
			     pcs->type == 42 ) && _nod_value_vector[k].find("SATURATION")
			    != string::npos)
				tec_file << ", SATURATION2";
		//-------------------------------------WW
	}
	const size_t mfp_value_vector_size (mfp_value_vector.size());
	for (size_t k = 0; k < mfp_value_vector_size; k++)
		//NB
		tec_file << ", \"" << mfp_value_vector[k] << "\"";

	// PCON
	//MX
	const size_t nPconName (_pcon_value_vector.size());
	for (size_t k = 0; k < nPconName; k++)
		//MX
		tec_file << ", " << _pcon_value_vector[k] << "";
	tec_file << endl;

	//--------------------------------------------------------------------
	// Write Header II: zone
	tec_file << "ZONE T=\"";
	tec_file << _time << "s\", ";
	//OK411
	tec_file << "N=" << m_msh->GetNodesNumber(false) << ", ";
	tec_file << "E=" << no_elements << ", ";
	tec_file << "F=" << "FEPOINT" << ", ";
	tec_file << "ET=" << e_type_name;
	tec_file << endl;
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   09/2004 OK Implementation
   01/2006 OK VAR,PCS,MSH concept
**************************************************************************/
void COutput::ELEWriteDOMDataTEC()
{
	//----------------------------------------------------------------------
	if(_ele_value_vector.empty ())
		return;
	//----------------------------------------------------------------------
	// File handling
	//......................................................................
	string tec_file_name = file_base_name + "_domain" + "_ele";
	if(getProcessType () != FiniteElement::INVALID_PROCESS) // PCS
		// 09/2010 TF msh_type_name;
		tec_file_name += "_" + convertProcessTypeToString (getProcessType());
	if(msh_type_name.size() > 1)          // MSH
		tec_file_name += "_" + msh_type_name;
	tec_file_name += TEC_FILE_EXTENSION;
	//WW
	if(!_new_file_opened)
		remove(tec_file_name.c_str());
	//......................................................................
	fstream tec_file (tec_file_name.data(),ios::app | ios::out);
	tec_file.setf(ios::scientific,ios::floatfield);
	tec_file.precision(12);
	if (!tec_file.good())
		return;
	tec_file.seekg(0L,ios::beg);
#ifdef SUPERCOMPUTER
	// kg44 buffer the output
	char mybuffer [MY_IO_BUFSIZE * MY_IO_BUFSIZE];
	tec_file.rdbuf()->pubsetbuf(mybuffer,MY_IO_BUFSIZE * MY_IO_BUFSIZE);
	//
#endif

	//--------------------------------------------------------------------
	WriteELEValuesTECHeader(tec_file);
	WriteELEValuesTECData(tec_file);
	//--------------------------------------------------------------------
	tec_file.close();                     // kg44 close file
}

void COutput::WriteELEValuesTECHeader(fstream &tec_file)
{
	// Write Header I: variables
	tec_file << "VARIABLES = \"X\",\"Y\",\"Z\",\"VX\",\"VY\",\"VZ\"";
	for (size_t i = 0; i < _ele_value_vector.size(); i++)
		//WW
		if (_ele_value_vector[i].find("VELOCITY") == string::npos)
			tec_file << "," << _ele_value_vector[i];
	tec_file << endl;

	// Write Header II: zone
	tec_file << "ZONE T=\"";
	tec_file << _time << "s\", ";
	tec_file << "I=" << (long) m_msh->ele_vector.size() << ", ";
	tec_file << "F=POINT" << ", ";
	tec_file << "C=BLACK";
	tec_file << endl;
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   09/2004 OK Implementation
   11/2005 OK MSH
   01/2006 OK
**************************************************************************/
void COutput::WriteELEValuesTECData(fstream &tec_file)
{
	CRFProcess* m_pcs = NULL;
	if (!_ele_value_vector.empty())
		m_pcs = GetPCS_ELE(_ele_value_vector[0]);
	else
		return;

	size_t no_ele_values = _ele_value_vector.size();
	bool out_element_vel = false;
	for (size_t j = 0; j < no_ele_values; j++) //WW

		if (_ele_value_vector[j].find("VELOCITY") != string::npos)
		{
			out_element_vel = true;
			break;
		}
	vector<int> ele_value_index_vector(no_ele_values);
	GetELEValuesIndexVector(ele_value_index_vector);

	MeshLib::CElem* m_ele = NULL;
	FiniteElement::ElementValue* gp_ele = NULL;
	for (size_t i = 0; i < m_msh->ele_vector.size(); i++)
	{
		m_ele = m_msh->ele_vector[i];
		double const* xyz(m_ele->GetGravityCenter());
		tec_file << xyz[0] << " " << xyz[1] << " " << xyz[2] << " ";
		if (out_element_vel)      //WW
		{
			if (PCSGet(FiniteElement::FLUID_MOMENTUM)) // PCH 16.11 2009
			{
				CRFProcess* pch_pcs = PCSGet(FiniteElement::FLUID_MOMENTUM);

				tec_file << pch_pcs->GetElementValue(i,
				                                     pch_pcs->GetElementValueIndex(
				                                             "VELOCITY1_X") + 1)
				<< " ";
				tec_file << pch_pcs->GetElementValue(i,
				                                     pch_pcs->GetElementValueIndex(
				                                             "VELOCITY1_Y") + 1)
				<< " ";
				tec_file << pch_pcs->GetElementValue(i,
				                                     pch_pcs->GetElementValueIndex(
				                                             "VELOCITY1_Z") + 1)
				<< " ";
			}
			else
			{
				gp_ele = ele_gp_value[i];
				tec_file << gp_ele->Velocity(0, 0) << " ";
				tec_file << gp_ele->Velocity(1, 0) << " ";
				tec_file << gp_ele->Velocity(2, 0) << " ";
			}
		}
		else
			for (size_t j = 0; j < ele_value_index_vector.size(); j++)
				tec_file
				<< m_pcs->GetElementValue(i, ele_value_index_vector[j])
				<< " ";
		/*
		   int j;
		   int eidx;
		   char ele_value_char[20];
		   int no_ele_values = (int)ele_value_vector.size();
		   for(j=0;j<no_ele_values;j++){
		   sprintf(ele_value_char,"%s",ele_value_vector[j].data());
		   eidx = PCSGetELEValueIndex(ele_value_char);
		   tec_file << ElGetElementVal(i,eidx) << " ";
		   }
		 */
		tec_file << endl;
	}

	ele_value_index_vector.clear();
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   08/2004 OK Implementation
   08/2004 WW Output by the order of distance to one end of the polyline
        OK this should be done in the PLY method
   08/2005 WW Changes due to the Geometry element class applied
           Remove existing files
   12/2005 OK Mass transport specifics
   12/2005 OK VAR,MSH,PCS concept
   12/2005 WW Output stress invariants
   08/2006 OK FLUX
   10/2008 OK MFP values
   07/2010 TF substituted GEOGetPLYByName
**************************************************************************/
double COutput::NODWritePLYDataTEC(int number)
{
	//WW  int nidx;
	long gnode;
	bool bdummy = false;
	int stress_i[6], strain_i[6];
	double ss[6];
	double val_n = 0.;                    //WW
	                                      //OK4704
	if ((_nod_value_vector.size() == 0) && (mfp_value_vector.size() == 0))
		return 0.0;

	// TF
	GEOLIB::Polyline const* const ply (
	        dynamic_cast<GEOLIB::Polyline const* const>(this->getGeoObj()));
	if (this->getGeoType() != GEOLIB::POLYLINE || ply == NULL)
	{
		std::cerr << "COutput::NODWritePLYDataTEC geometric object is not a polyline" <<
		std::endl;
		return 0.0;
	}

	// File handling
	std::string tec_file_name = file_base_name + "_ply_" + geo_name + "_t"
	                            + number2str<size_t> (_id); //OK4709
	if (getProcessType() != FiniteElement::INVALID_PROCESS)
		tec_file_name += "_" + convertProcessTypeToString(getProcessType());
	if (msh_type_name.size() > 0)
		tec_file_name += "_" + msh_type_name;
	tec_file_name += TEC_FILE_EXTENSION;
	if (!_new_file_opened)
		remove(tec_file_name.c_str());  //WW

	//WW
	fstream tec_file(tec_file_name.data(), ios::app | ios::out);

	tec_file.setf(ios::scientific, ios::floatfield);
	tec_file.precision(12);
	if (!tec_file.good())
		return 0.0;
	tec_file.seekg(0L, ios::beg);
#ifdef SUPERCOMPUTER
	// kg44 buffer the output
	char mybuffer [MY_IO_BUFSIZE * MY_IO_BUFSIZE];
	tec_file.rdbuf()->pubsetbuf(mybuffer,MY_IO_BUFSIZE * MY_IO_BUFSIZE);
	//
#endif
	//----------------------------------------------------------------------
	// Tests
	//......................................................................
	// GEO
//   CGLPolyline* m_ply = GEOGetPLYByName(geo_name);//CC
//   if (!m_ply)
//   {
//      cout << "Warning in COutput::NODWritePLYDataTEC - no GEO data" << endl;
//      tec_file << "Warning in COutput::NODWritePLYDataTEC - no GEO data: "
//         << geo_name << endl;
//      tec_file.close();
//      return 0.0;
//   }

	// MSH
	//	CFEMesh* m_msh = GetMSH();
	//	m_msh = GetMSH();
	if (!m_msh)
		cout << "Warning in COutput::NODWritePLYDataTEC - no MSH data" << endl;
	//OKtec_file << "Warning in COutput::NODWritePLYDataTEC - no MSH data: " << geo_name << endl;
	//OKtec_file.close();
	//OKToDo return;
	else
		m_msh->SwitchOnQuadraticNodes(false);  //WW

	// PCS
	if (getProcessType() == FiniteElement::INVALID_PROCESS)
		m_pcs = NULL;
	else
		m_pcs = PCSGet(getProcessType());

	CRFProcess* dm_pcs = NULL;            //WW
	for (size_t i = 0; i < pcs_vector.size(); i++)
		if (isDeformationProcess (pcs_vector[i]->getProcessType ()))
		{
			dm_pcs = pcs_vector[i];
			break;
		}

	/* //WW
	   // VEL
	   int v_eidx[3];
	   CRFProcess* m_pcs_flow (PCSGetFlow());
	   //m_pcs_flow = PCSGet("GROUNDWATER_FLOW"); //OKToDo
	   if (!m_pcs_flow)
	   {
	   //WW cout << "Warning in COutput::NODWritePLYDataTEC() - no PCS flow data" << endl;
	   //tec_file << "Warning in COutput::NODWritePLYDataTEC() - no PCS flow data " << endl;
	   //tec_file.close();
	   //return 0.0;
	   }
	   else
	   {
	   v_eidx[0] = m_pcs_flow->GetElementValueIndex("VELOCITY1_X");
	   v_eidx[1] = m_pcs_flow->GetElementValueIndex("VELOCITY1_Y");
	   v_eidx[2] = m_pcs_flow->GetElementValueIndex("VELOCITY1_Z");
	   }
	 */

//   for (size_t i = 0; i < 3; i++)
//   {
//      if (v_eidx[i] < 0)
//      {
//         //WW cout << "Warning in COutput::NODWritePLYDataTEC() - no PCS flow data" << endl;
//         //tec_file << "Warning in COutput::NODWritePLYDataTEC() - no PCS flow data " << endl;
//         //tec_file.close();
//      }
//   }
	//--------------------------------------------------------------------
	// NIDX for output variables
	size_t no_variables (_nod_value_vector.size());
	std::vector<int> NodeIndex(no_variables);
	GetNodeIndexVector(NodeIndex);
	//--------------------------------------------------------------------
	// Write header
	if (number == 0 || number == 1)       //WW if(number==1)
	{
		//project_title;
		std::string project_title_string = "Profiles along polylines";
		tec_file << " TITLE = \"" << project_title_string
		         << "\"" << endl;
		tec_file << " VARIABLES = \"DIST\" ";
		for (size_t k = 0; k < no_variables; k++)
		{
			tec_file << "\"" << _nod_value_vector[k] << "\" ";
			//-------------------------------------WW
			m_pcs = GetPCS(_nod_value_vector[k]);
			if (m_pcs && m_pcs->type == 1212 && _nod_value_vector[k].find(
			            "SATURATION") != string::npos)
				tec_file << "SATURATION2 ";
			//-------------------------------------WW
			if (_nod_value_vector[k].compare("FLUX") == 0)
				tec_file << "FLUX_INNER" << " ";
		}
		//....................................................................
		//OK4709
		for (size_t k = 0; k < mfp_value_vector.size(); k++)
			tec_file << "\"" << mfp_value_vector[k] << "\" ";
		//....................................................................
		// WW: M specific data
		if (dm_pcs)               //WW

			tec_file << " p_(1st_Invariant) " << " q_(2nd_Invariant)  "
			         << " Effective_Strain";
		tec_file << endl;
	}
	//....................................................................
	// WW: M specific data
	size_t ns = 4;
	if (dm_pcs)                           //WW
	{
		stress_i[0] = dm_pcs->GetNodeValueIndex("STRESS_XX");
		stress_i[1] = dm_pcs->GetNodeValueIndex("STRESS_YY");
		stress_i[2] = dm_pcs->GetNodeValueIndex("STRESS_ZZ");
		stress_i[3] = dm_pcs->GetNodeValueIndex("STRESS_XY");
		strain_i[0] = dm_pcs->GetNodeValueIndex("STRAIN_XX");
		strain_i[1] = dm_pcs->GetNodeValueIndex("STRAIN_YY");
		strain_i[2] = dm_pcs->GetNodeValueIndex("STRAIN_ZZ");
		strain_i[3] = dm_pcs->GetNodeValueIndex("STRAIN_XY");
		if (max_dim == 2)         // 3D
		{
			ns = 6;
			stress_i[4] = dm_pcs->GetNodeValueIndex("STRESS_XZ");
			stress_i[5] = dm_pcs->GetNodeValueIndex("STRESS_YZ");
			strain_i[4] = dm_pcs->GetNodeValueIndex("STRAIN_XZ");
			strain_i[5] = dm_pcs->GetNodeValueIndex("STRAIN_YZ");
		}
	}
	//......................................................................
	// , I=" << NodeListLength << ", J=1, K=1, F=POINT" << endl;
	tec_file << " ZONE T=\"TIME=" << _time << "\"" << endl;
	//----------------------------------------------------------------------
	// Write data
	//======================================================================
	double flux_sum = 0.0;                //OK
	double flux_nod;

	m_msh->SwitchOnQuadraticNodes(false); //WW
	// NOD at PLY
	std::vector<long> nodes_vector;

	CGLPolyline* m_ply = GEOGetPLYByName(geo_name);
//   m_msh->GetNODOnPLY(m_ply, old_nodes_vector); // TF

	double tmp_min_edge_length (m_msh->getMinEdgeLength());
	m_msh->setMinEdgeLength (m_ply->epsilon);
	m_msh->GetNODOnPLY(ply, nodes_vector); // TF
	std::vector<double> interpolation_points;
	m_msh->getPointsForInterpolationAlongPolyline (ply, interpolation_points);
	m_msh->setMinEdgeLength(tmp_min_edge_length);

//   std::cout << "size of nodes_vector: " << nodes_vector.size() << ", size of old_nodes_vector: " << old_nodes_vector.size() << std::endl;
	//bool b_specified_pcs = (m_pcs != NULL); //NW m_pcs = PCSGet(pcs_type_name);
	for (size_t j(0); j < nodes_vector.size(); j++)
	{
//		tec_file << m_ply->getSBuffer()[j] << " ";
		tec_file << interpolation_points[j] << " ";
		//WW
//		long old_gnode = nodes_vector[m_ply->getOrderedPoints()[j]];
		gnode = nodes_vector[j];
		for (size_t k = 0; k < no_variables; k++)
		{
			//if(!(_nod_value_vector[k].compare("FLUX")==0))  // removed JOD, does not work for multiple flow processes
			//if (!b_specified_pcs) //NW
			if (msh_type_name != "COMPARTMENT") // JOD 4.10.01
				m_pcs = PCSGet(_nod_value_vector[k], bdummy);

			if (!m_pcs)
			{
				cout << "Warning in COutput::NODWritePLYDataTEC - no PCS data"
				     << endl;
				tec_file
				<< "Warning in COutput::NODWritePLYDataTEC - no PCS data"
				<< endl;
				return 0.0;
			}
			// WW
//			double old_val_n = m_pcs->GetNodeValue(old_gnode, NodeIndex[k]);
			val_n = m_pcs->GetNodeValue(gnode, NodeIndex[k]);
//			tec_file << old_val_n << " ";
			tec_file << val_n << " ";
			if (m_pcs->type == 1212 && (_nod_value_vector[k].find("SATURATION")
			                            != string::npos))
				tec_file << 1. - val_n << " ";

			if (_nod_value_vector[k].compare("FLUX") == 0)
			{
				if (aktueller_zeitschritt == 0) //OK
					flux_nod = 0.0;
				else
					flux_nod = NODFlux(gnode);
				tec_file << flux_nod << " ";
				//flux_sum += abs(m_pcs->eqs->b[gnode]);
				flux_sum += abs(flux_nod);
				//OK cout << gnode << " " << flux_nod << " " << flux_sum << endl;
			}
		}
		if (dm_pcs) //WW
		{
			for (size_t i = 0; i < ns; i++)
				ss[i] = dm_pcs->GetNodeValue(gnode, stress_i[i]);
			tec_file << -DeviatoricStress(ss) / 3.0 << " ";
			tec_file << sqrt(3.0 * TensorMutiplication2(ss, ss,
			                                            m_msh->GetCoordinateFlag() /
			                                            10) / 2.0) << "  ";
			for (size_t i = 0; i < ns; i++)
				ss[i] = dm_pcs->GetNodeValue(gnode, strain_i[i]);
			DeviatoricStress(ss);
			tec_file << sqrt(3.0 * TensorMutiplication2(ss, ss,
			                                            m_msh->GetCoordinateFlag() /
			                                            10) / 2.0);
		}

		// MFP //OK4704
		//OK4704
		for (size_t k = 0; k < mfp_value_vector.size(); k++)
			//     tec_file << MFPGetNodeValue(gnode,mfp_value_vector[k],0) << " "; //NB
			tec_file << MFPGetNodeValue(gnode, mfp_value_vector[k], atoi(
			                                    &mfp_value_vector[k][mfp_value_vector[k
			                                                         ].size() - 1]) - 1)
			<< " ";  //NB: MFP output for all phases

		tec_file << endl;
	}
	tec_file.close();                     // kg44 close file
	return flux_sum;
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   08/2004 OK Implementation
   08/2005 WW MultiMesh
   12/2005 OK VAR,MSH,PCS concept
   12/2005 WW Output stress invariants
   10/2010 TF changed access to process type
**************************************************************************/
void COutput::NODWritePNTDataTEC(double time_current,int time_step_number)
{
	CRFProcess* dm_pcs = NULL;
	for (size_t i = 0; i < pcs_vector.size(); i++)
		//		if (pcs_vector[i]->pcs_type_name.find("DEFORMATION") != string::npos) { TF
		if (isDeformationProcess (pcs_vector[i]->getProcessType()))
		{
			dm_pcs = pcs_vector[i];
			break;
		}

	// File handling
	std::string tec_file_name(file_base_name + "_time_");
	addInfoToFileName(tec_file_name, true, true, true);

	if (!_new_file_opened)
		remove(tec_file_name.c_str());  //WW
	//......................................................................
	fstream tec_file(tec_file_name.data(), ios::app | ios::out);

	tec_file.setf(ios::scientific, ios::floatfield);
	tec_file.precision(12);
	if (!tec_file.good())
		return;
	tec_file.seekg(0L, ios::beg);
#ifdef SUPERCOMPUTER
	// kg44 buffer the output
	char mybuffer [MY_IO_BUFSIZE * MY_IO_BUFSIZE];
	tec_file.rdbuf()->pubsetbuf(mybuffer,MY_IO_BUFSIZE * MY_IO_BUFSIZE);
	//
#endif
	//--------------------------------------------------------------------
	// Tests
	//......................................................................
	//	CFEMesh* m_msh = GetMSH();
	//	m_msh = GetMSH();
	if (!m_msh)
	{
		cout << "Warning in COutput::NODWritePNTDataTEC - no MSH data: "
		     << endl;
		tec_file << "Warning in COutput::NODWritePNTDataTEC - no MSH data: "
		         << endl;
		tec_file.close();
		return;
	}

	//----------------------------------------------------------------------
	// NIDX for output variables
	size_t no_variables(_nod_value_vector.size());
	vector<int> NodeIndex(no_variables);
	GetNodeIndexVector(NodeIndex);
	//--------------------------------------------------------------------
	// Write header
	if (time_step_number == 0)            //WW  Old: if(time_step_number==1)
	{
		//project_title;
		std::string project_title_string = "Time curves in points";
		tec_file << " TITLE = \"" << project_title_string
		         << "\"" << endl;
		tec_file << " VARIABLES = \"TIME \" ";

		//    if(pcs_type_name.compare("RANDOM_WALK")==0)
		if (getProcessType() == FiniteElement::RANDOM_WALK)
			tec_file << "leavingParticles ";
		for (size_t k = 0; k < no_variables; k++) //WW
		{
			tec_file << " \"" << _nod_value_vector[k] << "\" ";
			//-------------------------------------WW
			m_pcs = GetPCS(_nod_value_vector[k]);
			if (m_pcs && m_pcs->type == 1212 && _nod_value_vector[k].find(
			            "SATURATION") != string::npos)
				tec_file << "SATURATION2 ";
			//-------------------------------------WW
		}
		//OK411
		for (size_t k = 0; k < mfp_value_vector.size(); k++)
			//NB MFP data names for multiple phases
			tec_file << " \"" << mfp_value_vector[k] << "\" ";
		//
#ifdef RFW_FRACTURE
		for(i = 0; i < (int)mmp_vector.size(); ++i)
			if( mmp_vector[i]->frac_num > 0)
				for(int j = 0; j < mmp_vector[i]->frac_num; ++j)
					tec_file << mmp_vector[i]->frac_names[j] << "_k " <<
					mmp_vector[i]->frac_names[j] << "_aper "
					         << mmp_vector[i]->frac_names[j] << "_closed ";

#endif

		if (dm_pcs)               //WW
			tec_file << " p_(1st_Invariant) " << " q_(2nd_Invariant)  "
			         << " Effective_Strain";
		tec_file << endl;

		if (geo_name.find("POINT") == std::string::npos)
			tec_file << " ZONE T=\"POINT="
			//, I=" << anz_zeitschritte << ", J=1, K=1, F=POINT" << endl;
			         << getGeoTypeAsString() << geo_name << "\"" << endl;
		else
			tec_file << " ZONE T=\"POINT=" << geo_name << "\""
			         << endl;  //, I=" << anz_zeitschritte << ", J=1, K=1, F=POINT" << endl;
	}

	// For deformation
	size_t ns = 4;
	int stress_i[6], strain_i[6];
	double ss[6];
	if (dm_pcs)                           //WW
	{
		stress_i[0] = dm_pcs->GetNodeValueIndex("STRESS_XX");
		stress_i[1] = dm_pcs->GetNodeValueIndex("STRESS_YY");
		stress_i[2] = dm_pcs->GetNodeValueIndex("STRESS_ZZ");
		stress_i[3] = dm_pcs->GetNodeValueIndex("STRESS_XY");
		strain_i[0] = dm_pcs->GetNodeValueIndex("STRAIN_XX");
		strain_i[1] = dm_pcs->GetNodeValueIndex("STRAIN_YY");
		strain_i[2] = dm_pcs->GetNodeValueIndex("STRAIN_ZZ");
		strain_i[3] = dm_pcs->GetNodeValueIndex("STRAIN_XY");
		if (m_msh->GetCoordinateFlag() / 10 == 3) // 3D
		{
			ns = 6;
			stress_i[4] = dm_pcs->GetNodeValueIndex("STRESS_XZ");
			stress_i[5] = dm_pcs->GetNodeValueIndex("STRESS_YZ");
			strain_i[4] = dm_pcs->GetNodeValueIndex("STRAIN_XZ");
			strain_i[5] = dm_pcs->GetNodeValueIndex("STRAIN_YZ");
		}
	}
	//--------------------------------------------------------------------
	// Write data
	//......................................................................
	tec_file << time_current << " ";
	//......................................................................
	// NOD values
	if (getProcessType() == FiniteElement::RANDOM_WALK)
		tec_file << m_msh->PT->leavingParticles << " ";
	int timelevel;
	CRFProcess* m_pcs_out = NULL;

	// fetch geometric entities, especial the associated GEOLIB::Point vector
	if (pcs_vector[0] == NULL)
		return;

	long msh_node_number(m_msh->GetNODOnPNT(
	                             static_cast<const GEOLIB::Point*> (getGeoObj())));

	// Mass transport
	if (getProcessType() == FiniteElement::MASS_TRANSPORT)
		for (size_t i = 0; i < _nod_value_vector.size(); i++)
		{
			std::string nod_value_name = _nod_value_vector[i];
			for (size_t l = 0; l < pcs_vector.size(); l++)
			{
				m_pcs = pcs_vector[l];
				//				if (m_pcs->pcs_type_name.compare("MASS_TRANSPORT") == 0) { TF
				if (m_pcs->getProcessType() == FiniteElement::MASS_TRANSPORT)
				{
					timelevel = 0;
					for (size_t m = 0; m < m_pcs->nod_val_name_vector.size();
					     m++)
						if (m_pcs->nod_val_name_vector[m].compare(
						            nod_value_name) == 0)
						{
							//							m_pcs_out = PCSGet(pcs_type_name, nod_value_name);
							m_pcs_out
							        = PCSGet(FiniteElement::MASS_TRANSPORT,
							                 nod_value_name);
							if (timelevel == 1)
							{
								int nidx =
								        m_pcs_out->
								        GetNodeValueIndex(
								                nod_value_name) +
								        timelevel;
								tec_file << m_pcs_out->GetNodeValue(
								        msh_node_number,
								        nidx) << " ";
							}
							timelevel++;
						}
				}
			}
		}
	else
	{
		double flux_nod, flux_sum = 0.0;
		for (size_t i = 0; i < _nod_value_vector.size(); i++)
		{
			// PCS
			if (!(_nod_value_vector[i].compare("FLUX") == 0)
				|| getProcessType() == FiniteElement::OVERLAND_FLOW) //JOD separate infiltration flux output in overland flow

				m_pcs = GetPCS(_nod_value_vector[i]);
			else
				m_pcs = GetPCS();
			if (!m_pcs)
			{
				cout << "Warning in COutput::NODWritePLYDataTEC - no PCS data"
				     << endl;
				tec_file
				<< "Warning in COutput::NODWritePLYDataTEC - no PCS data"
				<< endl;
				return;
			}
			//..................................................................
			// PCS
			if (!(_nod_value_vector[i].compare("FLUX") == 0)
				|| getProcessType() == FiniteElement::OVERLAND_FLOW) // JOD separate infiltration flux output in overland flow
			{
				//-----------------------------------------WW
				double val_n = m_pcs->GetNodeValue(msh_node_number,
				                                   NodeIndex[i]);
				tec_file << val_n << " ";
				m_pcs = GetPCS(_nod_value_vector[i]);
				if (m_pcs->type == 1212 && (_nod_value_vector[i].find(
				                                    "SATURATION") != string::npos))
					tec_file << 1. - val_n << " ";
				//-----------------------------------------WW
			}
			else
			{
				flux_nod = NODFlux(msh_node_number);
				tec_file << flux_nod << " ";
				//flux_sum += abs(m_pcs->eqs->b[gnode]);
				flux_sum += abs(flux_nod);
				//OK cout << gnode << " " << flux_nod << " " << flux_sum << endl;
			}
		}
		//....................................................................
#ifdef RFW_FRACTURE
		for(i = 0; i < (int)mmp_vector.size(); ++i)
			if( mmp_vector[i]->frac_num > 0)
				for(int j = 0; j < mmp_vector[i]->frac_num; ++j)
					tec_file << mmp_vector[i]->frac_perm[j] << " " <<
					mmp_vector[i]->avg_aperture[j] << " "
					         << mmp_vector[i]->closed_fraction[j] << " ";

#endif
		//....................................................................
		if (dm_pcs)               //WW
		{
			for (size_t i = 0; i < ns; i++)
				ss[i] = dm_pcs->GetNodeValue(msh_node_number, stress_i[i]);
			tec_file << -DeviatoricStress(ss) / 3.0 << " ";
			tec_file << sqrt(3.0 * TensorMutiplication2(ss, ss,
			                                            m_msh->GetCoordinateFlag() /
			                                            10) / 2.0) << "  ";
			for (size_t i = 0; i < ns; i++)
				ss[i] = dm_pcs->GetNodeValue(msh_node_number, strain_i[i]);
			DeviatoricStress(ss);
			tec_file << sqrt(3.0 * TensorMutiplication2(ss, ss,
			                                            m_msh->GetCoordinateFlag() /
			                                            10) / 2.0);
		}
		//OK411
		for (size_t k = 0; k < mfp_value_vector.size(); k++)
			tec_file << MFPGetNodeValue(msh_node_number, mfp_value_vector[k],
			                            atoi(&mfp_value_vector[k][mfp_value_vector[k].
			                                                      size() - 1])
			                            - 1) << " ";  //NB
	}
	tec_file << endl;
	//----------------------------------------------------------------------
	tec_file.close();                     // kg44 close file
}

void COutput::WriteRFOHeader(fstream &rfo_file)
{
	//#0#0#0#1#0.00000000000000e+000#0#3915###########################################
	rfo_file << "#0#0#0#1#";
	rfo_file << _time;
	rfo_file << "#0#";
	rfo_file << OGS_VERSION;
	rfo_file << "###########################################";
	rfo_file << endl;
}

void COutput::WriteRFONodes(fstream &rfo_file)
{
	//0 101 100
	rfo_file << 0 << " " << m_msh->nod_vector.size() << " "
	         << m_msh->ele_vector.size() << std::endl;
	//0 0.00000000000000 0.00000000000000 0.00000000000000
	for (size_t i = 0; i < m_msh->nod_vector.size(); i++)
	{
		double const* const pnt_i (m_msh->nod_vector[i]->getData());
		rfo_file << i << " " << pnt_i[0] << " " << pnt_i[1] << " "
		         << pnt_i[2] << " " << std::endl;
	}
}

void COutput::WriteRFOElements(fstream &rfo_file)
{
	size_t j;
	MeshLib::CElem* m_ele = NULL;
	//0 0 -1 line 0 1
	for(long i = 0; i < (long)m_msh->ele_vector.size(); i++)
	{
		m_ele = m_msh->ele_vector[i];
		rfo_file << i << " " << m_ele->GetPatchIndex() \
		         << " -1 " \
		         << m_ele->GetName() << " ";
		for(j = 0; j < m_ele->GetNodesNumber(false); j++)
			rfo_file << m_ele->getNodeIndices()[j] << " ";
		rfo_file << endl;
	}
}

void COutput::WriteRFOValues(fstream &rfo_file)
{
	int p,nidx;
	CRFProcess* m_pcs = NULL;
	//1 2 4
	rfo_file << "1 1 4" << endl;
	//2 1 1
	int no_processes = (int)pcs_vector.size();
	rfo_file << no_processes;
	for(p = 0; p < no_processes; p++)
		rfo_file << " 1";
	rfo_file << endl;
	//PRESSURE1, Pa
	// Names and units
	for(p = 0; p < no_processes; p++)
	{
		m_pcs = pcs_vector[p];
		rfo_file << m_pcs->pcs_primary_function_name[0];
		rfo_file << ", ";
		rfo_file << m_pcs->pcs_primary_function_unit[0];
		rfo_file << endl;
	}
	//0 0.00000000000000 0.00000000000000
	// Node values
	for(long i = 0; i < (long)m_msh->nod_vector.size(); i++)
	{
		rfo_file << i;
		for(p = 0; p < no_processes; p++)
		{
			m_pcs = pcs_vector[p];
			nidx = m_pcs->GetNodeValueIndex(m_pcs->pcs_primary_function_name[0]) + 1;
			rfo_file << " " << m_pcs->GetNodeValue(i,nidx);
		}
		rfo_file << " " << endl;
	}
}

void COutput::WriteRFO()
{
	m_msh = FEMGet(convertProcessTypeToString (getProcessType()));
	if(!m_msh)
	{
		cout << "Warning in COutput::WriteRFONodes - no MSH data" << endl;
		return;
	}
	//--------------------------------------------------------------------
	// File handling
	string rfo_file_name;
	rfo_file_name = file_base_name + "." + "rfo";
	//WW
	if(!_new_file_opened)
		remove(rfo_file_name.c_str());
	fstream rfo_file (rfo_file_name.data(),ios::app | ios::out);

	rfo_file.setf(ios::scientific,ios::floatfield);
	rfo_file.precision(12);
	if (!rfo_file.good())
		return;
	rfo_file.seekg(0L,ios::beg);
#ifdef SUPERCOMPUTER
	// kg44 buffer the output
	char mybuffer [MY_IO_BUFSIZE * MY_IO_BUFSIZE];
	rfo_file.rdbuf()->pubsetbuf(mybuffer,MY_IO_BUFSIZE * MY_IO_BUFSIZE);
	//
#endif
	//--------------------------------------------------------------------
	WriteRFOHeader(rfo_file);
	WriteRFONodes(rfo_file);
	WriteRFOElements(rfo_file);
	WriteRFOValues(rfo_file);
	//  RFOWriteELEValues();
	rfo_file.close();                     // kg44 close file
}

void COutput::NODWriteSFCDataTEC(int number)
{
	if (_nod_value_vector.size() == 0)
		return;

	m_msh = FEMGet(convertProcessTypeToString(getProcessType()));
	m_pcs = PCSGet(getProcessType());

	// File handling
	char number_char[3];
	sprintf(number_char, "%i", number);
	string number_string = number_char;
	//	string tec_file_name = pcs_type_name + "_sfc_" + geo_name + "_t"
	//				+ number_string + TEC_FILE_EXTENSION;
	std::string tec_file_name = convertProcessTypeToString (getProcessType()) + "_sfc_" +
	                            geo_name + "_t"
	                            + number_string + TEC_FILE_EXTENSION;
	if (!_new_file_opened)
		remove(tec_file_name.c_str());  //WW
	fstream tec_file(tec_file_name.data(), ios::app | ios::out);
	tec_file.setf(ios::scientific, ios::floatfield);
	tec_file.precision(12);
	if (!tec_file.good())
		return;
	tec_file.seekg(0L, ios::beg);
#ifdef SUPERCOMPUTER
	// kg44 buffer the output
	char mybuffer [MY_IO_BUFSIZE * MY_IO_BUFSIZE];
	tec_file.rdbuf()->pubsetbuf(mybuffer,MY_IO_BUFSIZE * MY_IO_BUFSIZE);
	//
#endif
	//--------------------------------------------------------------------
	// Write header
	//project_title;
	string project_title_string = "Profile at surface";
	tec_file << " TITLE = \"" << project_title_string << "\""
	         << endl;
	tec_file << " VARIABLES = \"X\",\"Y\",\"Z\",";
	for (size_t k = 0; k < _nod_value_vector.size(); k++)
		tec_file << _nod_value_vector[k] << ",";
	tec_file << endl;
	// , I=" << NodeListLength << ", J=1, K=1, F=POINT" << endl;
	tec_file << " ZONE T=\"TIME=" << _time << "\"" << endl;
	//--------------------------------------------------------------------
	// Write data
	std::vector<long> nodes_vector;
	Surface* m_sfc = NULL;
	m_sfc = GEOGetSFCByName(geo_name);    //CC
	if (m_sfc)
	{
		m_msh->GetNODOnSFC(m_sfc, nodes_vector);
		for (size_t i = 0; i < m_msh->nod_vector.size(); i++)
		{
			double const* const pnt_i (m_msh->nod_vector[i]->getData());
			tec_file << pnt_i[0] << " ";
			tec_file << pnt_i[1] << " ";
			tec_file << pnt_i[2] << " ";
			for (size_t k = 0; k < _nod_value_vector.size(); k++)
			{
				int nidx = m_pcs->GetNodeValueIndex(_nod_value_vector[k]) + 1;
				tec_file << m_pcs->GetNodeValue(nodes_vector[i], nidx) << " ";
			}
			tec_file << endl;
		}
	}
	else
		tec_file << "Error in NODWritePLYDataTEC: polyline " << geo_name
		         << " not found" << endl;
	tec_file.close();                     // kg44 close file
}

/**************************************************************************
   FEMLib-Method:
   12/2005 OK Implementation
   04/2006 CMCD no mesh option & flux weighting
   last modification:
**************************************************************************/
void COutput::NODWriteSFCAverageDataTEC(double time_current,int time_step_number)
{
	bool no_pcs = false;
	double dtemp;
	vector<long> sfc_nodes_vector;
	double node_flux = 0.0;
	int idx = -1;
	double t_flux = 0.0;
	double node_conc = 0.0;
	CRFProcess* m_pcs_gw (PCSGet(FiniteElement::GROUNDWATER_FLOW));
	if (!m_pcs_gw)
		PCSGet(FiniteElement::LIQUID_FLOW);
	//--------------------------------------------------------------------
	// Tests
	Surface* m_sfc = NULL;
	m_sfc = GEOGetSFCByName(geo_name);
	if (!m_sfc)
	{
		cout << "Warning in COutput::NODWriteSFCAverageDataTEC - no GEO data"
		     << endl;
		return;
	}
	//	CFEMesh* m_msh = NULL;
	m_msh = FEMGet(convertProcessTypeToString(getProcessType()));
	if (!m_msh)
	{
		cout << "Warning in COutput::NODWriteSFCAverageDataTEC - no MSH data"
		     << endl;
		return;
	}
	CRFProcess* m_pcs(PCSGet(getProcessType()));
	if (!m_pcs)
		no_pcs = true;
	//cout << "Warning in COutput::NODWriteSFCAverageDataTEC - no PCS data" << endl;
	//return;
	//--------------------------------------------------------------------
	// File handling
	string tec_file_name = file_base_name + "_TBC_" + getGeoTypeAsString()
	                       + "_" + geo_name + TEC_FILE_EXTENSION;
	if (!_new_file_opened)
		remove(tec_file_name.c_str());  //WW
	fstream tec_file(tec_file_name.data(), ios::app | ios::out);
	tec_file.setf(ios::scientific, ios::floatfield);
	tec_file.precision(12);
	if (!tec_file.good())
		return;
	tec_file.seekg(0L, ios::beg);
#ifdef SUPERCOMPUTER
	// kg44 buffer the output
	char mybuffer [MY_IO_BUFSIZE * MY_IO_BUFSIZE];
	tec_file.rdbuf()->pubsetbuf(mybuffer,MY_IO_BUFSIZE * MY_IO_BUFSIZE);
	//
#endif
	//--------------------------------------------------------------------
	// Write header
	if (time_step_number == 0)            // WW Old:  if(time_step_number==1)
	{
		//project_title;
		string project_title_string = "Time curve at surface";
		tec_file << " TITLE = \"" << project_title_string
		         << "\"" << endl;
		tec_file << " VARIABLES = Time ";
		for (size_t i = 0; i < _nod_value_vector.size(); i++)
			tec_file << _nod_value_vector[i] << " ";
		tec_file << endl;
		//, I=" << anz_zeitschritte << ", J=1, K=1, F=POINT" << endl;
		tec_file << " ZONE T=\"SFC=" << geo_name << "\"" << endl;
	}
	//--------------------------------------------------------------------
	// node_value_index_vector

	std::vector<int> node_value_index_vector(_nod_value_vector.size());
	//	int itemp;
	if (m_pcs)
		for (size_t i = 0; i < _nod_value_vector.size(); i++)
		{
			node_value_index_vector[i] = m_pcs->GetNodeValueIndex(
			        _nod_value_vector[i]) + 1;
			//			itemp = node_value_index_vector[i];
			for (size_t n_pv = 0; n_pv < m_pcs->GetPrimaryVNumber(); n_pv++)
				if (_nod_value_vector[i].compare(
				            m_pcs->pcs_primary_function_name[n_pv]) == 0)
				{
					node_value_index_vector[i]++;
					break;
				}
		}
	//--------------------------------------------------------------------
	// Write data
	if (no_pcs)
	{
		tec_file << time_current << " ";
		for (size_t i = 0; i < _nod_value_vector.size(); i++)
		{
			//Specified currently for MASS_TRANSPORT only.
			m_pcs = PCSGet(FiniteElement::MASS_TRANSPORT, _nod_value_vector[i]);
			node_value_index_vector[i] = m_pcs->GetNodeValueIndex(
			        _nod_value_vector[i]) + 1;
			m_pcs->m_msh->GetNODOnSFC(m_sfc, sfc_nodes_vector);
			dtemp = 0.0;
			t_flux = 0.0;
			for (size_t j = 0; j < sfc_nodes_vector.size(); j++)
			{
				idx = m_pcs_gw->GetNodeValueIndex("FLUX");
				node_flux = abs(
				        m_pcs_gw->GetNodeValue(sfc_nodes_vector[j], idx));
				node_conc = m_pcs->GetNodeValue(sfc_nodes_vector[j],
				                                node_value_index_vector[i]);
				dtemp += (node_flux * node_conc);
				t_flux += node_flux;
			}
			dtemp /= t_flux;
			tec_file << dtemp << " ";
		}
		tec_file << endl;
	}
	else
	{
		m_msh->GetNODOnSFC(m_sfc, sfc_nodes_vector);
		idx = m_pcs_gw->GetNodeValueIndex("FLUX");
		tec_file << time_current << " ";
		dtemp = 0.0;
		t_flux = 0.0;
		for (size_t i = 0; i < _nod_value_vector.size(); i++)
		{
			dtemp = 0.0;
			for (size_t j = 0; j < sfc_nodes_vector.size(); j++)
			{
				node_flux = abs(
				        m_pcs_gw->GetNodeValue(sfc_nodes_vector[j], idx));
				node_conc = m_pcs->GetNodeValue(sfc_nodes_vector[j],
				                                node_value_index_vector[i]);
				dtemp += (node_flux * node_conc);
				t_flux += node_flux;
			}
			dtemp /= t_flux;
			tec_file << dtemp << " ";
		}
		tec_file << endl;
	}
	tec_file.close();                     // kg44 close file
}

void COutput::GetNodeIndexVector(vector<int>&NodeIndex)
{
	CRFProcess* pcs = NULL;
	const size_t nName = _nod_value_vector.size();
	if (getProcessType() != FiniteElement::INVALID_PROCESS)
	{
		pcs = PCSGet(getProcessType());
		for (size_t k = 0; k < nName; k++)
		{
			if (getProcessType () == FiniteElement::MASS_TRANSPORT)
				pcs = PCSGet(getProcessType(), _nod_value_vector[k]);
			if (!pcs)
			{
				cout
				<< "Warning in COutput::GetNodeIndexVector - no PCS data: "
				<< _nod_value_vector[k] << endl;
				return;
			}
			NodeIndex[k] = pcs->GetNodeValueIndex(_nod_value_vector[k],true);  // JT latest
		}
	}
	else if (msh_type_name.size() > 0)
	{
		pcs = PCSGet(msh_type_name);
		if (!pcs)
		{
			cout << "Warning in COutput::GetNodeIndexVector - no PCS data"
			     << endl;
			return;
		}
		for (size_t k = 0; k < nName; k++)
		{
			NodeIndex[k] = pcs->GetNodeValueIndex(_nod_value_vector[k],true); // JT latest
		}
	}
	else if (fem_msh_vector.size() == 1)
	{
		bool bdummy = true;
		for (size_t k = 0; k < nName; k++)
		{
			pcs = PCSGet(_nod_value_vector[k], bdummy);
			if (!pcs)
			{
				cout
				<< "Warning in COutput::GetNodeIndexVector - no PCS data: "
				<< _nod_value_vector[k] << endl;
				return;
			}
			NodeIndex[k] = pcs->GetNodeValueIndex(_nod_value_vector[k],true); // JT latest
		}
	}
}

CRFProcess* COutput::GetPCS(const string &var_name)
{
	CRFProcess* m_pcs = NULL;
	if(getProcessType () != FiniteElement::INVALID_PROCESS)
		m_pcs = PCSGet(getProcessType());
	else if(msh_type_name.size() > 0)
		m_pcs = PCSGet(msh_type_name);
	if(!m_pcs)
		m_pcs = PCSGet(var_name,true);
	return m_pcs;
}

// 09/2010 TF
CRFProcess* COutput::GetPCS()
{
	if(getProcessType () != FiniteElement::INVALID_PROCESS)
	{
		if (getProcess() != NULL)
			return getProcess();
		else
			return PCSGet(getProcessType());
	}
	else
	{
		CRFProcess* pcs (NULL);
		if(msh_type_name.size() > 0)
			pcs = PCSGet(msh_type_name);
		//	  if(!pcs)
		//		pcs = PCSGet(var_name,true);
		return pcs;
	}
}

CRFProcess* COutput::GetPCS_ELE(const string &var_name)
{
	string pcs_var_name;
	CRFProcess* m_pcs = NULL;
	//----------------------------------------------------------------------
	//  if(pcs_type_name.size()>0)
	//    m_pcs = PCSGet(pcs_type_name);
	if (getProcessType() != FiniteElement::INVALID_PROCESS)
		m_pcs = PCSGet(getProcessType());
	else if (msh_type_name.size() > 0)
		m_pcs = PCSGet(msh_type_name);

	if (!m_pcs)
		for (size_t i = 0; i < pcs_vector.size(); i++)
		{
			m_pcs = pcs_vector[i];
			for (int j = 0; j < m_pcs->pcs_number_of_evals; j++)
			{
				pcs_var_name = m_pcs->pcs_eval_name[j];
				if (pcs_var_name.compare(var_name) == 0)
					return m_pcs;
			}
		}
	return m_pcs;
}

void COutput::GetELEValuesIndexVector(vector<int>&ele_value_index_vector)
{
	if (_ele_value_vector[0].size() == 0)
		return;
	CRFProcess* m_pcs(GetPCS_ELE(_ele_value_vector[0]));
	for (size_t i = 0; i < _ele_value_vector.size(); i++)
		ele_value_index_vector[i] = m_pcs->GetElementValueIndex(_ele_value_vector[i]);
}

/**************************************************************************
   FEMLib-Method:
   Programing:
   01/2006 OK Implementation
   07/2010 TF substituted GEOGetPLYByName, renamed local variables for CGLPolyline,
         CFEMesh and CRFProcess
         (conflicting with class attributes that have the same name)
**************************************************************************/
void COutput::SetNODFluxAtPLY()
{
	// Tests
	//  CGLPolyline* ply = GEOGetPLYByName(geo_name);
	//	CGLPolyline* ply = polyline_vector[getGeoObjIdx()];
	//	if (!ply) {
	//		cout << "Warning in COutput::SetNODFluxAtPLY() - no PLY data" << endl;
	//		return;
	//	}

	//	CFEMesh* msh = GetMSH();
	if (!m_msh)
	{
		cout << "Warning in COutput::SetNODFluxAtPLY() - no MSH data" << endl;
		return;
	}

	CRFProcess* pcs = GetPCS("FLUX");
	if (!pcs)
	{
		cout << "Warning in COutput::SetNODFluxAtPLY() - no PCS data" << endl;
		return;
	}

	std::vector<long> nodes_vector;
	//	msh->GetNODOnPLY(ply, nodes_vector);
	m_msh->GetNODOnPLY(static_cast<const GEOLIB::Polyline*>(getGeoObj()), nodes_vector);
	std::vector<double> node_value_vector;
	node_value_vector.resize(nodes_vector.size());
	int nidx = pcs->GetNodeValueIndex("FLUX");
	for (size_t i = 0; i < nodes_vector.size(); i++)
		node_value_vector[i] = pcs->GetNodeValue(nodes_vector[i], nidx);
	//----------------------------------------------------------------------
	//m_st->FaceIntegration(m_pcs,nodes_vector,node_value_vector);
}

void COutput::ELEWriteSFC_TEC()
{
	//----------------------------------------------------------------------
	if (_ele_value_vector.size() == 0)
		return;
	//----------------------------------------------------------------------
	// File handling
	//......................................................................
	std::string tec_file_name = file_base_name + "_surface" + "_ele";
	addInfoToFileName(tec_file_name, false, true, true);
	//  if(pcs_type_name.size()>1) // PCS
	//    tec_file_name += "_" + msh_type_name;
	//  if(msh_type_name.size()>1) // MSH
	//    tec_file_name += "_" + msh_type_name;
	//  tec_file_name += TEC_FILE_EXTENSION;

	if (!_new_file_opened)
		remove(tec_file_name.c_str());  //WW
	//......................................................................
	fstream tec_file(tec_file_name.data(), ios::app | ios::out);
	tec_file.setf(ios::scientific, ios::floatfield);
	tec_file.precision(12);
	if (!tec_file.good())
		return;
	tec_file.seekg(0L, ios::beg);
#ifdef SUPERCOMPUTER
	// kg44 buffer the output
	char mybuffer [MY_IO_BUFSIZE * MY_IO_BUFSIZE];
	tec_file.rdbuf()->pubsetbuf(mybuffer,MY_IO_BUFSIZE * MY_IO_BUFSIZE);
	//
#endif
	//--------------------------------------------------------------------
	vector<long> tmp_ele_sfc_vector;
	tmp_ele_sfc_vector.clear();
	//--------------------------------------------------------------------
	ELEWriteSFC_TECHeader(tec_file);
	ELEWriteSFC_TECData(tec_file);
	//--------------------------------------------------------------------
	tec_file.close();                     // kg44 close file
}

void COutput::ELEWriteSFC_TECHeader(fstream &tec_file)
{
	// Write Header I: variables
	tec_file << "VARIABLES = \"X\",\"Y\",\"Z\"";
	for (size_t i = 0; i < _ele_value_vector.size(); i++)
		tec_file << "," << _ele_value_vector[i];
	tec_file << endl;
	//--------------------------------------------------------------------
	// Write Header II: zone
	tec_file << "ZONE T=\"";
	tec_file << _time << "s\", ";
	tec_file << "I=" << (long) m_msh->ele_vector.size() << ", ";
	tec_file << "F=POINT" << ", ";
	tec_file << "C=BLACK";
	tec_file << endl;
}

void COutput::ELEWriteSFC_TECData(fstream &tec_file)
{
	tec_file << "COutput::ELEWriteSFC_TECData - implementation not finished" << endl;

	/* // Make it as comment to avoid compilation warnings. 18.082011 WW
	   long i;
	   int j;
	   MeshLib::CElem* m_ele = NULL;
	   MeshLib::CElem* m_ele_neighbor = NULL;
	   double v[3];
	   CRFProcess* m_pcs = NULL;
	   double v_n;
	   //--------------------------------------------------------------------
	   m_pcs = pcs_vector[0];                         //GetPCS_ELE(ele_value_vector[0]);
	   int nidx[3];
	   nidx[0] = m_pcs->GetElementValueIndex("VELOCITY1_X");
	   nidx[1] = m_pcs->GetElementValueIndex("VELOCITY1_Y");
	   nidx[2] = m_pcs->GetElementValueIndex("VELOCITY1_Z");
	   //--------------------------------------------------------------------
	   for(i=0l;i<(long)m_msh->ele_vector.size();i++)
	   {
	   m_ele = m_msh->ele_vector[i];
	   for(j=0;j<m_ele->GetFacesNumber();j++)
	   {
	      m_ele_neighbor = m_ele->GetNeighbor(j);
	      if((m_ele->GetDimension() - m_ele_neighbor->GetDimension())==1)
	      {
	         v[0] = m_pcs->GetElementValue(m_ele->GetIndex(),nidx[0]);
	         v[1] = m_pcs->GetElementValue(m_ele->GetIndex(),nidx[1]);
	         v[2] = m_pcs->GetElementValue(m_ele->GetIndex(),nidx[2]);
	         m_ele_neighbor->SetNormalVector();

	         v_n = v[0]*m_ele_neighbor->normal_vector[0] \
	 + v[1]*m_ele_neighbor->normal_vector[1] \
	 + v[2]*m_ele_neighbor->normal_vector[2];

	      }
	   }
	   }
	 */
	//--------------------------------------------------------------------
}

/**************************************************************************
   FEMLib-Method:
   08/2006 OK Implementation
   modification:
   05/2010 TF if geo_type_name[3] == 'N' (case point) only a pointer to the
   CGLPoint is fetched, but the pointer is never used
   07/2010 TF substituted GEOGetPLYByName
   09/2010 TF substituted pcs_type_name
**************************************************************************/
void COutput::CalcELEFluxes()
{
	const FiniteElement::ProcessType pcs_type (getProcessType());
	if (pcs_type == FiniteElement::INVALID_PROCESS)      // WW moved it here.

		//WW cout << "Warning in COutput::CalcELEFluxes(): no PCS data" << endl;
		return;

	CRFProcess* pcs = PCSGet(getProcessType());
	if (FiniteElement::isDeformationProcess(pcs_type) || !isFlowProcess (pcs_type)
	    //WW
	    || pcs->m_msh->geo_name.find("REGIONAL") != string::npos)
		return;
	//----------------------------------------------------------------------
	switch (getGeoType())
	{
	case GEOLIB::POLYLINE:
	{
		//		CGLPolyline* ply = GEOGetPLYByName(geo_name);
		//		if (!ply)
		//			std::cout << "Warning in COutput::CalcELEFluxes - no GEO data" << std::endl;
		double f_n_sum = 0.0;
		//		f_n_sum = pcs->CalcELEFluxes(ply); // TF
		f_n_sum = pcs->CalcELEFluxes(static_cast<const GEOLIB::Polyline*> (getGeoObj()));

		ELEWritePLY_TEC();
		//BUGFIX_4402_OK_1
		TIMValue_TEC(f_n_sum);
		break;
	}
	case GEOLIB::SURFACE:
	{
		//		Surface* m_sfc = GEOGetSFCByName(geo_name);
		//		pcs->CalcELEFluxes(m_sfc);
		break;
	}
	case GEOLIB::VOLUME:
	{
		//		CGLVolume* m_vol = GEOGetVOL(geo_name);
		//		pcs->CalcELEFluxes(m_vol);
		break;
	}
	case GEOLIB::GEODOMAIN:               //domAin
		//pcs->CalcELEFluxes(m_dom);
		break;
	default:
		cout << "Warning in COutput::CalcELEFluxes(): no GEO type data" << endl;
	}

	// WW   pcs->CalcELEFluxes(ply) changed 'mark' of elements
	for (size_t i = 0; i < fem_msh_vector.size(); i++)
		for (size_t j = 0; j < fem_msh_vector[i]->ele_vector.size(); j++)
			fem_msh_vector[i]->ele_vector[j]->MarkingAll(true);
}

void COutput::ELEWritePLY_TEC()
{
	//----------------------------------------------------------------------
	if(_ele_value_vector.size() == 0)
		return;
	//----------------------------------------------------------------------
	// File handling
	//......................................................................
	string tec_file_name = file_base_name; // + "_ply" + "_ele";
	tec_file_name += "_" + getGeoTypeAsString();
	tec_file_name += "_" + geo_name;
	tec_file_name += "_ELE";
	//  if(pcs_type_name.size()>1) // PCS
	//    tec_file_name += "_" + pcs_type_name;
	if(getProcessType () != FiniteElement::INVALID_PROCESS) // PCS
		tec_file_name += "_" + convertProcessTypeToString (getProcessType ());

	if(msh_type_name.size() > 1)          // MSH
		tec_file_name += "_" + msh_type_name;
	tec_file_name += TEC_FILE_EXTENSION;
	if(!_new_file_opened)
		remove(tec_file_name.c_str());  //WW
	//......................................................................
	fstream tec_file(tec_file_name.data(),ios::app | ios::out);
	tec_file.setf(ios::scientific,ios::floatfield);
	tec_file.precision(12);
	if (!tec_file.good())
		return;
	tec_file.seekg(0L,ios::beg);
#ifdef SUPERCOMPUTER
	// kg44 buffer the output
	char mybuffer [MY_IO_BUFSIZE * MY_IO_BUFSIZE];
	tec_file.rdbuf()->pubsetbuf(mybuffer,MY_IO_BUFSIZE * MY_IO_BUFSIZE);
	//
#endif
	//--------------------------------------------------------------------
	vector<long>tmp_ele_ply_vector;
	tmp_ele_ply_vector.clear();
	//--------------------------------------------------------------------
	ELEWritePLY_TECHeader(tec_file);
	ELEWritePLY_TECData(tec_file);
	//--------------------------------------------------------------------

	tec_file.close();                     // kg44 close file
}

void COutput::ELEWritePLY_TECHeader(fstream &tec_file)
{
	// Write Header I: variables
	tec_file << "VARIABLES = \"X\",\"Y\",\"Z\"";
	for (size_t i = 0; i < _ele_value_vector.size(); i++)
		tec_file << "," << _ele_value_vector[i];
	tec_file << endl;

	// Write Header II: zone
	tec_file << "ZONE T=\"";
	tec_file << _time << "s\", ";
	tec_file << endl;
}

/**************************************************************************
   FEMLib-Method:
   06/2006 OK Implementation
   07/2010 TF substituted GEOGetPLYByName
   10/2010 TF changed access to process type
**************************************************************************/
void COutput::ELEWritePLY_TECData(fstream &tec_file)
{
	//	CRFProcess* pcs = PCSGet(pcs_type_name);
	CRFProcess* pcs = PCSGet(getProcessType());
	int f_eidx[3];
	f_eidx[0] = pcs->GetElementValueIndex("FLUX_X");
	f_eidx[1] = pcs->GetElementValueIndex("FLUX_Y");
	f_eidx[2] = pcs->GetElementValueIndex("FLUX_Z");
	for (size_t i = 0; i < 3; i++)
		if (f_eidx[i] < 0)
		{
			cout <<
			"Fatal error in CRFProcess::CalcELEFluxes(CGLPolyline*m_ply) - abort";
			abort();
		}
	int v_eidx[3];
	CRFProcess* m_pcs_flow = NULL;

	//	if (pcs->pcs_type_name.find("FLOW") != string::npos) {
	if (isFlowProcess(pcs->getProcessType ()))
		m_pcs_flow = pcs;
	else
		m_pcs_flow = PCSGet(FiniteElement::GROUNDWATER_FLOW);

	v_eidx[0] = m_pcs_flow->GetElementValueIndex("VELOCITY1_X");
	v_eidx[1] = m_pcs_flow->GetElementValueIndex("VELOCITY1_Y");
	v_eidx[2] = m_pcs_flow->GetElementValueIndex("VELOCITY1_Z");
	for (size_t i = 0; i < 3; i++)
		if (v_eidx[i] < 0)
		{
			cout <<
			"Fatal error in CRFProcess::CalcELEFluxes(CGLPolyline*m_ply) - abort";
			abort();
		}

	//  CGLPolyline* m_ply = GEOGetPLYByName(geo_name);

	// Get elements at GEO
	//	vector<long> ele_vector_at_geo;
	//	m_msh->GetELEOnPLY(m_ply, ele_vector_at_geo);
	std::vector<size_t> ele_vector_at_geo;
	m_msh->GetELEOnPLY(static_cast<const GEOLIB::Polyline*> (getGeoObj()), ele_vector_at_geo);

	// helper variables
	Math_Group::vec<MeshLib::CEdge*> ele_edges_vector(15);
	Math_Group::vec<MeshLib::CNode*> edge_nodes(3);
	double edge_mid_vector[3] = {0.0,0.0,0.0};

	for (size_t i = 0; i < ele_vector_at_geo.size(); i++)
	{
		MeshLib::CElem* m_ele = m_msh->ele_vector[ele_vector_at_geo[i]];
		// x,y,z
		m_ele->GetEdges(ele_edges_vector);
		for (size_t j = 0; j < m_ele->GetEdgesNumber(); j++)
		{
			MeshLib::CEdge* m_edg = ele_edges_vector[j];
			if (m_edg->GetMark())
			{
				m_edg->GetNodes(edge_nodes);
				double const* const pnt0(edge_nodes[0]->getData());
				double const* const pnt1(edge_nodes[1]->getData());
				edge_mid_vector[0] = 0.5 * (pnt1[0] + pnt0[0]);
				edge_mid_vector[1] = 0.5 * (pnt1[1] + pnt0[1]);
				edge_mid_vector[2] = 0.5 * (pnt1[2] + pnt0[2]);
			}
		}
		tec_file << edge_mid_vector[0] << " " << edge_mid_vector[1] << " "
		         << edge_mid_vector[2];
		// ele vector values
		tec_file << " " << m_pcs_flow->GetElementValue(m_ele->GetIndex(),
		                                               v_eidx[0]);
		tec_file << " " << m_pcs_flow->GetElementValue(m_ele->GetIndex(),
		                                               v_eidx[1]);
		tec_file << " " << m_pcs_flow->GetElementValue(m_ele->GetIndex(),
		                                               v_eidx[2]);
		tec_file << " " << pcs->GetElementValue(m_ele->GetIndex(), f_eidx[0]);
		tec_file << " " << pcs->GetElementValue(m_ele->GetIndex(), f_eidx[1]);
		tec_file << " " << pcs->GetElementValue(m_ele->GetIndex(), f_eidx[2]);
		tec_file << endl;
	}
}

void COutput::TIMValue_TEC(double tim_value)
{
	//----------------------------------------------------------------------
	// File handling
	//......................................................................
	fstream tec_file;
	string tec_file_name = file_base_name; // + "_ply" + "_ele";
	tec_file_name += "_" + getGeoTypeAsString();
	tec_file_name += "_" + geo_name;
	tec_file_name += "_TIM";
	//  if(pcs_type_name.size()>1) // PCS
	//    tec_file_name += "_" + pcs_type_name;
	if(getProcessType () != FiniteElement::INVALID_PROCESS) // PCS
		tec_file_name += "_" + convertProcessTypeToString (getProcessType());
	if(msh_type_name.size() > 1)          // MSH
		tec_file_name += "_" + msh_type_name;
	tec_file_name += TEC_FILE_EXTENSION;
	if(!_new_file_opened)
		remove(tec_file_name.c_str());  //WW
	//......................................................................
	tec_file.open(tec_file_name.data(),ios::app | ios::out);
	tec_file.setf(ios::scientific,ios::floatfield);
	tec_file.precision(12);
	if (!tec_file.good())
		return;
	tec_file.seekg(0L,ios::beg);
#ifdef SUPERCOMPUTER
	// kg44 buffer the output
	char mybuffer [MY_IO_BUFSIZE * MY_IO_BUFSIZE];
	tec_file.rdbuf()->pubsetbuf(mybuffer,MY_IO_BUFSIZE * MY_IO_BUFSIZE);
	//
#endif
	//--------------------------------------------------------------------
	// Write Header I: variables
	if(aktueller_zeitschritt == 1)
	{
		tec_file << "VARIABLES = \"Time\",\"Value\"";
		tec_file << endl;
		//--------------------------------------------------------------------
		// Write Header II: zone
		tec_file << "ZONE T=";
		tec_file << geo_name;
		tec_file << endl;
	}
	//--------------------------------------------------------------------
	tec_file << aktuelle_zeit << " " << tim_value << endl;
	//--------------------------------------------------------------------
	tec_file.close();                     // kg44 close file
}

double COutput::NODFlux(long nod_number)
{
	nod_number = nod_number;              //OK411
	/*
	   cout << gnode << " " \
	    << m_pcs->GetNodeValue(gnode,NodeIndex[k]) << end
	   flux_sum += m_pcs->GetNodeValue(gnode,NodeIndex[k]);
	 */
	// All elements at node //OK
#ifdef NEW_EQS                                 //WW. 07.11.2008
	return 0.;                            //To do: m_pcs->eqs_new->b[nod_number];
#else
	// Element nodal RHS contributions
	m_pcs->eqs->b[nod_number] = 0.0;
	MeshLib::CNode* nod = m_msh->nod_vector[nod_number];
	m_pcs->AssembleParabolicEquationRHSVector(nod);
	return m_pcs->eqs->b[nod_number];
#endif
}

void COutput::NODWriteLAYDataTEC(int time_step_number)
{
	// Tests
	const size_t nName (_nod_value_vector.size());
	if (nName == 0)
		return;
	std::vector<int> NodeIndex(nName);

	// PCS
	CRFProcess* m_pcs = PCSGet(getProcessType());
	if (!m_pcs)
		return;
	for (size_t k = 0; k < nName; k++)
		NodeIndex[k] = m_pcs->GetNodeValueIndex(_nod_value_vector[k]);

	// MSH
	//	m_msh = GetMSH();
	if (!m_msh)
	{
		cout << "Warning in COutput::NODWriteLAYDataTEC() - no MSH data"
		     << endl;
		return;
	}

	// File name handling
	char char_time_step_number[10];
	sprintf(char_time_step_number, "%i", time_step_number);
	string tec_file_name(file_base_name);
	tec_file_name += "_layer_";
	tec_file_name += char_time_step_number;
	tec_file_name += TEC_FILE_EXTENSION;
	fstream tec_file(tec_file_name.data(), ios::trunc | ios::out);

	tec_file.setf(ios::scientific, ios::floatfield);
	tec_file.precision(12);
	if (!tec_file.good())
		return;
#ifdef SUPERCOMPUTER
	//
	// kg44 buffer the output
	char mybuffer [MY_IO_BUFSIZE * MY_IO_BUFSIZE];
	tec_file.rdbuf()->pubsetbuf(mybuffer,MY_IO_BUFSIZE * MY_IO_BUFSIZE);
	//
#endif
	//--------------------------------------------------------------------
	// Write Header I: variables
	tec_file << "VARIABLES = X,Y,Z,N";
	for (size_t k = 0; k < nName; k++)
		tec_file << "," << _nod_value_vector[k] << " ";
	tec_file << endl;

	long j;
	long no_per_layer = m_msh->GetNodesNumber(false)
	                    / (m_msh->getNumberOfMeshLayers() + 1);
	long jl;
	for (size_t l = 0; l < m_msh->getNumberOfMeshLayers() + 1; l++)
	{
		//--------------------------------------------------------------------
		tec_file << "ZONE T=LAYER" << l << endl;
		//--------------------------------------------------------------------
		for (j = 0l; j < no_per_layer; j++)
		{
			jl = j + j* m_msh->getNumberOfMeshLayers() + l;
			//..................................................................
			// XYZ
			double const* const pnt (m_msh->nod_vector[jl]->getData());
			tec_file << pnt[0] << " ";
			tec_file << pnt[1] << " ";
			tec_file << pnt[2] << " ";
			tec_file << jl << " ";
			//..................................................................
			for (size_t k = 0; k < nName; k++)
				tec_file << m_pcs->GetNodeValue(
				        m_msh->nod_vector[jl]->GetIndex(), NodeIndex[k]) << " ";
			tec_file << endl;
		}
	}
	tec_file.close();                     // kg44 close file
}

/**************************************************************************
   FEMLib-Method:
   Task: Write PCON data for ChemApp output
   Programing:
   08/2008 MX Implementation
**************************************************************************/
void COutput::PCONWriteDOMDataTEC()
{
	int te = 0;
	string eleType;
	string tec_file_name;
#if defined(USE_MPI) || defined(USE_MPI_PARPROC) || defined(USE_MPI_REGSOIL)
	char tf_name[10];
	std::cout << "Process " << myrank << " in WriteDOMDataTEC" << "\n";
#endif

	//----------------------------------------------------------------------
	// Tests
	if(_pcon_value_vector.size() == 0)
		return;
	//......................................................................
	// MSH
	//m_msh = FEMGet(pcs_type_name);
	//  m_msh = GetMSH();
	if(!m_msh)
	{
		cout << "Warning in COutput::NODWriteDOMDataTEC() - no MSH data" << endl;
		return;
	}
	//======================================================================
	vector<int> mesh_type_list;           //NW
	if (m_msh->getNumberOfLines() > 0)
		mesh_type_list.push_back(1);
	if (m_msh->getNumberOfQuads () > 0)
		mesh_type_list.push_back(2);
	if (m_msh->getNumberOfHexs () > 0)
		mesh_type_list.push_back(3);
	if (m_msh->getNumberOfTris () > 0)
		mesh_type_list.push_back(4);
	if (m_msh->getNumberOfTets () > 0)
		mesh_type_list.push_back(5);
	if (m_msh->getNumberOfPrisms () > 0)
		mesh_type_list.push_back(6);

	// Output files for each mesh type
	//NW
	for (int i = 0; i < (int)mesh_type_list.size(); i++)
	{
		te = mesh_type_list[i];
		//----------------------------------------------------------------------
		// File name handling
		tec_file_name = file_base_name + "_" + "domain_PCON";
		if(msh_type_name.size() > 0) // MultiMSH
			tec_file_name += "_" + msh_type_name;
		//  if(pcs_type_name.size()>0) // PCS
		//    tec_file_name += "_" + pcs_type_name;
		if(getProcessType () != FiniteElement::INVALID_PROCESS) // PCS
			tec_file_name += "_" + convertProcessTypeToString (getProcessType());
		//======================================================================
		switch (te)               //NW
		{
		case 1:
			tec_file_name += "_line";
			eleType = "QUADRILATERAL";
			break;
		case 2:
			tec_file_name += "_quad";
			eleType = "QUADRILATERAL";
			break;
		case 3:
			tec_file_name += "_hex";
			eleType = "BRICK";
			break;
		case 4:
			tec_file_name += "_tri";
			eleType = "QUADRILATERAL";
			break;
		case 5:
			tec_file_name += "_tet";
			eleType = "TETRAHEDRON";
			break;
		case 6:
			tec_file_name += "_pris";
			eleType = "BRICK";
			break;
		}
		/*
		   if(m_msh->msh_no_line>0)
		   {
		      tec_file_name += "_line";
		      eleType = "QUADRILATERAL";
		     te=1;
		   }
		   else if (m_msh->msh_no_quad>0)
		   {
		      tec_file_name += "_quad";
		      eleType = "QUADRILATERAL";
		   te=2;
		   }
		   else if (m_msh->msh_no_hexs>0)
		   {
		   tec_file_name += "_hex";
		   eleType = "BRICK";
		   te=3;
		   }
		   else if (m_msh->msh_no_tris>0)
		   {
		   tec_file_name += "_tri";
		   //???Who was this eleType = "TRIANGLE";
		   eleType = "QUADRILATERAL";
		   te=4;
		   }
		   else if (m_msh->msh_no_tets>0)
		   {
		   tec_file_name += "_tet";
		   eleType = "TETRAHEDRON";
		   te=5;
		   }
		   else if (m_msh->msh_no_pris>0)
		   {
		   tec_file_name += "_pris";
		   eleType = "BRICK";
		   te=6;
		   }
		 */
#if defined(USE_MPI) || defined(USE_MPI_PARPROC) || defined(USE_MPI_REGSOIL)
		sprintf(tf_name, "%d", myrank);
		tec_file_name += "_" + string(tf_name);
		std::cout << "Tecplot filename: " << tec_file_name << endl;
#endif
		tec_file_name += TEC_FILE_EXTENSION;
		//WW
		if(!_new_file_opened)
			remove(tec_file_name.c_str());
		fstream tec_file (tec_file_name.data(),ios::app | ios::out);
		tec_file.setf(ios::scientific,ios::floatfield);
		tec_file.precision(12);
		if (!tec_file.good())
			return;
#ifdef SUPERCOMPUTER
		// kg44 buffer the output
		char mybuf1 [MY_IO_BUFSIZE * MY_IO_BUFSIZE];
		tec_file.rdbuf()->pubsetbuf(mybuf1,MY_IO_BUFSIZE * MY_IO_BUFSIZE);
#endif
		//
		WriteTECHeader(tec_file,te,eleType);
		WriteTECNodePCONData(tec_file);
		WriteTECElementData(tec_file,te);
		tec_file.close();         // kg44 close file
		//--------------------------------------------------------------------
		// tri elements
		// ***** 07/2010 TF commented out block since the global variable is always zero
		//    if(msh_no_tris>0){
		//    //string tec_file_name = pcs_type_name + "_" + "domain" + "_tri" + TEC_FILE_EXTENSION;
		//#ifdef SUPERCOMPUTER
		//// buffer the output
		//      char sxbuf1[MY_IO_BUFSIZE*MY_IO_BUFSIZE];
		//#endif
		//      string tec_file_name = file_base_name + "_" + "domain" + "_tri" + TEC_FILE_EXTENSION;
		//      fstream tec_file1 (tec_file_name.data(),ios::app|ios::out);
		//      tec_file1.setf(ios::scientific,ios::floatfield);
		//      tec_file1.precision(12);
		//      if (!tec_file1.good()) return;
		//#ifdef SUPERCOMPUTER
		//      tec_file1.rdbuf()->pubsetbuf(sxbuf1,MY_IO_BUFSIZE*MY_IO_BUFSIZE);
		//#endif
		//      //OK  tec_file1.clear();
		//      //OK  tec_file1.seekg(0L,ios::beg);
		//      WriteTECHeader(tec_file1,4,"TRIANGLE");
		//      WriteTECNodeData(tec_file1);
		//      WriteTECElementData(tec_file1,4);
		//      tec_file1.close(); // kg44 close file
		//    }
		//--------------------------------------------------------------------
		// ***** 07/2010 TF commented out block since the global variable is always zero
		// quad elements
		//    if(msh_no_quad>0){
		//      //string tec_file_name = pcs_type_name + "_" + "domain" + "_quad" + TEC_FILE_EXTENSION;
		//#ifdef SUPERCOMPUTER
		//      char sxbuf2[MY_IO_BUFSIZE*MY_IO_BUFSIZE];
		//#endif
		//      string tec_file_name = file_base_name + "_" + "domain" + "_quad" + TEC_FILE_EXTENSION;
		//      fstream tec_file (tec_file_name.data(),ios::app|ios::out);
		//
		//      tec_file.setf(ios::scientific,ios::floatfield);
		//      tec_file.precision(12);
		//      if (!tec_file.good()) return;
		//#ifdef SUPERCOMPUTER
		//      tec_file.rdbuf()->pubsetbuf(sxbuf2,MY_IO_BUFSIZE*MY_IO_BUFSIZE);
		//#endif
		//      WriteTECHeader(tec_file,2,"QUADRILATERAL");
		//      WriteTECNodeData(tec_file);
		//      WriteTECElementData(tec_file,2);
		//      tec_file.close(); // kg44 close file
		//    }
		//--------------------------------------------------------------------
		// ***** 07/2010 TF commented out block since the global variable is always zero
		// tet elements
		//    if(msh_no_tets>0){
		//      //string tec_file_name = pcs_type_name + "_" + "domain" + "_tet" + TEC_FILE_EXTENSION;
		//#ifdef SUPERCOMPUTER
		//      char sxbuf3[MY_IO_BUFSIZE*MY_IO_BUFSIZE];
		//#endif
		//
		//      string tec_file_name = file_base_name + "_" + "domain" + "_tet";
		//
		//#if defined(USE_MPI) || defined(USE_MPI_PARPROC) || defined(USE_MPI_REGSOIL)
		//      sprintf(tf_name, "%d", myrank);
		//      tec_file_name += "_" + string(tf_name);
		//#endif
		//
		//      tec_file_name += TEC_FILE_EXTENSION;
		//
		//      fstream tec_file (tec_file_name.data(),ios::app|ios::out);
		//
		//      tec_file.setf(ios::scientific,ios::floatfield);
		//      tec_file.precision(12);
		//      if (!tec_file.good()) return;
		//#ifdef SUPERCOMPUTER
		//      tec_file.rdbuf()->pubsetbuf(sxbuf3,MY_IO_BUFSIZE*MY_IO_BUFSIZE);
		//#endif
		//
		//      WriteTECHeader(tec_file,5,"TETRAHEDRON");
		//      WriteTECNodeData(tec_file);
		//      WriteTECElementData(tec_file,5);
		//      tec_file.close(); // kg44 close file
		//    }
		//--------------------------------------------------------------------
		// ***** 07/2010 TF commented out block since the global variable is always zero
		// pris elements
		//    if(msh_no_pris>0){
		//      //string tec_file_name = pcs_type_name + "_" + "domain" + "_pris" + TEC_FILE_EXTENSION;
		//#ifdef SUPERCOMPUTER
		//        char sxbuf4[MY_IO_BUFSIZE*MY_IO_BUFSIZE];
		//#endif
		//      string tec_file_name = file_base_name + "_" + "domain" + "_pris" + TEC_FILE_EXTENSION;
		//      fstream tec_file (tec_file_name.data(),ios::app|ios::out);
		//
		//      tec_file.setf(ios::scientific,ios::floatfield);
		//      tec_file.precision(12);
		//      if (!tec_file.good()) return;
		//#ifdef SUPERCOMPUTER
		//      tec_file.rdbuf()->pubsetbuf(sxbuf4,MY_IO_BUFSIZE*MY_IO_BUFSIZE);
		//#endif
		//
		//      WriteTECHeader(tec_file,6,"BRICK");
		//      WriteTECNodeData(tec_file);
		//      WriteTECElementData(tec_file,6);
		//      tec_file.close(); // kg44 close file
		//    }
		//--------------------------------------------------------------------
		// ***** 07/2010 TF commented out block since the global variable is always zero
		// hex elements
		//    if(msh_no_hexs>0){
		//      //string tec_file_name = pcs_type_name + "_" + "domain" + "_hex" + TEC_FILE_EXTENSION;
		//#ifdef SUPERCOMPUTER
		//        char sxbuf5[MY_IO_BUFSIZE*MY_IO_BUFSIZE];
		//#endif
		//
		//      string tec_file_name = file_base_name + "_" + "domain" + "_hex" + TEC_FILE_EXTENSION;
		//      fstream tec_file (tec_file_name.data(),ios::app|ios::out);
		//
		//
		//      tec_file.setf(ios::scientific,ios::floatfield);
		//      tec_file.precision(12);
		//      if (!tec_file.good()) return;
		//#ifdef SUPERCOMPUTER
		//      tec_file.rdbuf()->pubsetbuf(sxbuf5,MY_IO_BUFSIZE*MY_IO_BUFSIZE);
		//#endif
		//      WriteTECHeader(tec_file,3,"BRICK");
		//      WriteTECNodeData(tec_file);
		//      WriteTECElementData(tec_file,3);
		//      tec_file.close(); // kg44 close file
		//    }
	}
}

/**************************************************************************
   FEMLib-Method:
   Task: Node value output of PCON in aquous
   Programing:
   08/2008 MX Implementation
**************************************************************************/
void COutput::WriteTECNodePCONData(fstream &tec_file)
{
	const size_t nName (_pcon_value_vector.size());
	int nidx_dm[3];
	std::vector<int> PconIndex(nName);

	//  m_msh = GetMSH();

#ifdef CHEMAPP
	CEqlink* eq = NULL;

	eq = eq->GetREACTION();
	if (!eq)
		return;
	const int nPCON_aq = eq->NPCON[1];    //GetNPCON(1);
	eq->GetPconNameAq();

	for(i = 0; i < nName; i++)
	{
		for(k = 0; k < nPCON_aq; k++)
			//		 pcon_value_name = PconName_Aq[i];
			if(pcon_value_vector[i].compare(PconName_Aq[k]) == 0)
			{
				PconIndex[i] = k;
				break;
			}
	}
#endif
	// MSH
	//--------------------------------------------------------------------
	for (size_t j = 0l; j < m_msh->GetNodesNumber(false); j++)
	{
		// XYZ
		double x[3] =
		{m_msh->nod_vector[j]->getData()[0], m_msh->nod_vector[j]->getData()[1],
		 m_msh->nod_vector[j]->getData()[2]};
//      x[0] = m_msh->nod_vector[j]->X();
//      x[1] = m_msh->nod_vector[j]->Y();
//      x[2] = m_msh->nod_vector[j]->Z();
		// Amplifying DISPLACEMENTs
		if (M_Process || MH_Process) //WW

			for (size_t k = 0; k < max_dim + 1; k++)
				x[k] += out_amplifier * m_pcs->GetNodeValue(
				        m_msh->nod_vector[j]->GetIndex(), nidx_dm[k]);
		for (size_t i = 0; i < 3; i++)
			tec_file << x[i] << " ";
		// NOD values
#ifdef CHEMAPP
		for(size_t k = 0; k < nName; k++)
			tec_file << eq->GetPconAq_mol_amount(j,PconIndex[k]) << " ";

#endif
		tec_file << endl;
	}
}

void COutput::checkConsistency ()
{
	if (!_nod_value_vector.empty())
	{
		std::vector<std::string> del_index, alias_del_lindex;
		bool found = false;
		CRFProcess* pcs = NULL;
		const size_t pcs_vector_size(pcs_vector.size());
		const size_t nod_value_vector_size (_nod_value_vector.size());
		for (size_t j = 0; j < nod_value_vector_size; j++)
		{
			found = false; // initialize variable found
			for (size_t l = 0; l < pcs_vector_size; l++)
			{
				pcs = pcs_vector[l];
				for (size_t m = 0; m < pcs->nod_val_name_vector.size(); m++)
					if (pcs->nod_val_name_vector[m].compare(
					            _nod_value_vector[j]) == 0)
					{
						found = true;
						del_index.push_back(_nod_value_vector[j]);
                        alias_del_lindex.push_back(_alias_nod_value_vector[j]);
						break;
					}
				// end for(m...)
			}             // end for(l...)
			if (!found)
			{
				std::cout << "Warning - no PCS data for output variable "
				          << _nod_value_vector[j] << " in ";
				switch (getGeoType())
				{
				case GEOLIB::POINT:
					std::cout << "POINT " << getGeoName() << std::endl;
					break;
				case GEOLIB::POLYLINE:
					std::cout << "POLYLINE " << getGeoName() << std::endl;
					break;
				case GEOLIB::SURFACE:
					std::cout << "SURFACE " << getGeoName() << std::endl;
					break;
				case GEOLIB::VOLUME:
					std::cout << "VOLUME " << getGeoName() << std::endl;
					break;
				case GEOLIB::GEODOMAIN:
					std::cout << "DOMAIN " << getGeoName() << std::endl;
					break;
				case GEOLIB::INVALID:
					std::cout <<
					"WARNING: COutput::checkConsistency - invalid geo type" <<
					endl;
					break;
				}
			}
		}                         // end for(j...)

		// Reduce vector out->_nod_value_vector by elements which have no PCS
		if (del_index.size() < _nod_value_vector.size())
		{
			std::cout << " Reducing output to variables with existing PCS-data " <<
			std::endl;
			_nod_value_vector.clear();
			for (size_t j = 0; j < del_index.size(); j++)
				_nod_value_vector.push_back(del_index[j]);
            _alias_nod_value_vector.clear();
            for (size_t j = 0; j < del_index.size(); j++)
                _alias_nod_value_vector.push_back(alias_del_lindex[j]);
		}
		if (!pcs)
			pcs = this->GetPCS();
		if (!pcs)
			cout << "Warning in OUTData - no PCS data" << endl;
	}                                     // end if(_nod_value_vector.size()>0)
}

/**************************************************************************
   FEMLib-Method:
   Task: Set output variable names for internal use
   Programing:
   11/2011 NW Implementation
**************************************************************************/
void COutput::setInternalVarialbeNames(CFEMesh *msh)
{
#if 0
    if (_alias_nod_value_vector.empty())
        return;
    bool isXZplane = (msh->GetCoordinateFlag()==22);
    bool isPVD = (dat_type_name.compare("PVD") == 0); //currently only for PVD

    if (isXZplane && isPVD) {
        std::cout << "-> recognized XZ plane for PVD output." << std::endl;
        map<string,string> map_output_variable_name;
        map_output_variable_name.insert(pair<string, string>("DISPLACEMENT_Y1", "DISPLACEMENT_Z1" ));
        map_output_variable_name.insert(pair<string, string>("DISPLACEMENT_Z1", "DISPLACEMENT_Y1" ));
        map_output_variable_name.insert(pair<string, string>("STRESS_XY", "STRESS_XZ" ));
        map_output_variable_name.insert(pair<string, string>("STRESS_YY", "STRESS_ZZ" ));
        map_output_variable_name.insert(pair<string, string>("STRESS_ZZ", "STRESS_YY" ));
        map_output_variable_name.insert(pair<string, string>("STRESS_XZ", "STRESS_XY" ));
        map_output_variable_name.insert(pair<string, string>("STRAIN_XY", "STRAIN_XZ" ));
        map_output_variable_name.insert(pair<string, string>("STRAIN_YY", "STRAIN_ZZ" ));
        map_output_variable_name.insert(pair<string, string>("STRAIN_ZZ", "STRAIN_YY" ));
        map_output_variable_name.insert(pair<string, string>("STRAIN_XZ", "STRAIN_XY"  ));
        map_output_variable_name.insert(pair<string, string>("VELOCITY_Y1", "VELOCITY_Z1"));
        map_output_variable_name.insert(pair<string, string>("VELOCITY_Z1", "VELOCITY_Y1"));
        map_output_variable_name.insert(pair<string, string>("VELOCITY_Y2", "VELOCITY_Z2"));
        map_output_variable_name.insert(pair<string, string>("VELOCITY_Z2", "VELOCITY_Y2"));

        for (size_t j = 0; j < _alias_nod_value_vector.size(); j++) {
            if (map_output_variable_name.count(_alias_nod_value_vector[j])>0) {
                _nod_value_vector.push_back(map_output_variable_name[_alias_nod_value_vector[j]]);
            } else {
                _nod_value_vector.push_back(_alias_nod_value_vector[j]);
            }
        }
    } else {
        for (size_t j = 0; j < _alias_nod_value_vector.size(); j++) {
            _nod_value_vector.push_back(_alias_nod_value_vector[j]);
        }
    }
#else
    if (_nod_value_vector.empty())
        return;
    bool isXZplane = (msh->GetCoordinateFlag()==22);
    bool isPVD = (dat_type_name.compare("PVD") == 0); //currently only for PVD

    if (isXZplane && isPVD) {
        std::cout << "-> recognized XZ plane for PVD output." << std::endl;
        map<string,string> map_output_variable_name;
        map_output_variable_name.insert(pair<string, string>("DISPLACEMENT_Y1", "DISPLACEMENT_Z1" ));
        map_output_variable_name.insert(pair<string, string>("DISPLACEMENT_Z1", "DISPLACEMENT_Y1" ));
        map_output_variable_name.insert(pair<string, string>("STRESS_XY", "STRESS_XZ" ));
        map_output_variable_name.insert(pair<string, string>("STRESS_YY", "STRESS_ZZ" ));
        map_output_variable_name.insert(pair<string, string>("STRESS_ZZ", "STRESS_YY" ));
        map_output_variable_name.insert(pair<string, string>("STRESS_XZ", "STRESS_XY" ));
        map_output_variable_name.insert(pair<string, string>("STRAIN_XY", "STRAIN_XZ" ));
        map_output_variable_name.insert(pair<string, string>("STRAIN_YY", "STRAIN_ZZ" ));
        map_output_variable_name.insert(pair<string, string>("STRAIN_ZZ", "STRAIN_YY" ));
        map_output_variable_name.insert(pair<string, string>("STRAIN_XZ", "STRAIN_XY"  ));
        map_output_variable_name.insert(pair<string, string>("VELOCITY_Y1", "VELOCITY_Z1"));
        map_output_variable_name.insert(pair<string, string>("VELOCITY_Z1", "VELOCITY_Y1"));
        map_output_variable_name.insert(pair<string, string>("VELOCITY_Y2", "VELOCITY_Z2"));
        map_output_variable_name.insert(pair<string, string>("VELOCITY_Z2", "VELOCITY_Y2"));

        for (size_t j = 0; j < _nod_value_vector.size(); j++) {
            if (map_output_variable_name.count(_nod_value_vector[j])>0) {
                _alias_nod_value_vector.push_back(map_output_variable_name[_nod_value_vector[j]]);
            } else {
                _alias_nod_value_vector.push_back(_nod_value_vector[j]);
            }
        }
    } else {
        for (size_t j = 0; j < _nod_value_vector.size(); j++) {
            _alias_nod_value_vector.push_back(_nod_value_vector[j]);
        }
    }
#endif
}

void COutput::addInfoToFileName (std::string& file_name, bool geo, bool process, bool mesh) const
{
	// add geo type name
	if (geo)
		//		file_name += getGeoTypeAsString();
		file_name += geo_name;

	// add process type name
	if (getProcessType() != FiniteElement::INVALID_PROCESS && process)
		file_name += "_" + FiniteElement::convertProcessTypeToString (getProcessType());

	// add mesh type name
	if (msh_type_name.size() > 0 && mesh)
		file_name += "_" + msh_type_name;

	// finally add file extension
	file_name += TEC_FILE_EXTENSION;
}
