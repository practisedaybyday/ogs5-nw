/**************************************************************************
   FEMLib - Object: OUT
   Task:
   Programing:
   06/2004 OK Implementation
   last modified:
**************************************************************************/
#include "Configure.h"
#include "makros.h"
// C++ STL
#include <cfloat>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <list>
#include <string>
using namespace std;

// Base
#include "StringTools.h"

// FEM-Makros
#include "files0.h"
#include "makros.h"
// GeoSys-GeoLib
#include "GEOObjects.h"
#include "files0.h"
#include "geo_ply.h"
#include "geo_sfc.h"
// GeoSys-FEMLib
#include "LegacyVtkInterface.h"
#include "Output.h"
#include "fem_ele_std.h"
#include "mathlib.h"
#include "rf_msp_new.h"
#include "rf_out_new.h"
#include "rf_pcs.h"
#include "rf_pcs.h"
#include "rf_random_walk.h"
#include "rf_tim_new.h"
// GeoSys-MSHLib
#include "msh_lib.h"

// FileIO/FEMIO
#include "FEMIO/GeoIO.h"

#include "problem.h"

// Base
#include "StringTools.h"

extern size_t max_dim;                            //OK411 todo

#ifdef CHEMAPP
#include "eqlink.h"
#endif
#include "vtk.h"
// MPI Parallel
#if defined(USE_MPI) || defined(USE_MPI_PARPROC) || defined(USE_MPI_REGSOIL)
#include "par_ddc.h"
#endif
#ifdef SUPERCOMPUTER
// kg44 this is usefull for io-buffering as endl flushes the buffer
#define endl '\n'
#define MY_IO_BUFSIZE 4096
#endif
#ifdef GEM_REACT
#include "rf_REACT_GEM.h"
#endif
using MeshLib::CFEMesh;
//==========================================================================
vector<COutput*>out_vector;

/**************************************************************************
   FEMLib-Method:
   Task: OUT read function
   Programing:
   06/2004 OK Implementation
   08/2004 WW Remove the old files
   01/2005 OK Boolean type
   01/2005 OK Destruct before read
   06/2006 WW Remove the old files by new way
   06/2010 TF reformated, restructured, signature changed, use new GEOLIB data structures
**************************************************************************/
bool OUTRead(const std::string& file_base_name,
             const GEOLIB::GEOObjects& geo_obj,
             const std::string& unique_name)
{
	char line[MAX_ZEILE];
	std::string line_string;
	ios::pos_type position;
	bool output_version = false; // 02.2011. WW

	// File handling
	std::string out_file_name = file_base_name + OUT_FILE_EXTENSION;
	std::ifstream out_file(out_file_name.data(), ios::in);
	if (!out_file.good())
		return false;
	out_file.seekg(0L, ios::beg);

	// Keyword loop
	cout << "OUTRead" << endl;
	while (!out_file.eof())
	{
		out_file.getline(line, MAX_ZEILE);
		line_string = line;
		if (line_string.find("#STOP") != string::npos)
			return true;

		COutput* out(new COutput(out_vector.size()));
		out->getFileBaseName() = file_base_name;
		// Give version in file name
		//15.01.2008. WW
		if (line_string.find("#VERSION") != string::npos)
			output_version = true;  // 02.2011. WW
		//----------------------------------------------------------------------
		// keyword found
		if (line_string.find("#OUTPUT") != string::npos)
		{
			position = out->Read(out_file, geo_obj, unique_name);

			if(output_version) //// 02.2011. WW
			{
				std::string VersionStr = OGS_VERSION; //02.2011 WX
				int curPos = 0;
				int pos = 0;
				while((pos = VersionStr.find("/",curPos)) != -1)
				{
					VersionStr.replace(pos, 1, "_");
					curPos = pos + 1;
				}
				out->getFileBaseName().append("(V");
				out->getFileBaseName().append(VersionStr);
				out->getFileBaseName().append(")");
			}

			out_vector.push_back(out);

			//			char number_char[3]; //OK4709
			//			sprintf(number_char, "%i", (int) out_vector.size() - 1); //OK4709
			//			out->ID = number_char; //OK4709
			//			out->setID (out_vector.size() - 1);

			out_file.seekg(position, ios::beg);
		}                         // keyword found
	}                                     // eof
	return true;
}

/**************************************************************************
   FEMLib-Method: OUTWrite
   Task: master write function
   Programing:
   06/2004 OK Implementation
   last modification:
**************************************************************************/
void OUTWrite(string base_file_name)
{
	//========================================================================
	// File handling
	string out_file_name = base_file_name + OUT_FILE_EXTENSION;
	fstream out_file (out_file_name.data(),ios::trunc | ios::out);
	out_file.setf(ios::scientific,ios::floatfield);
	out_file.precision(12);
	if (!out_file.good())
		return;
	out_file.seekg(0L,ios::beg);
#ifdef SUPERCOMPUTER
	// kg44 buffer the output
	char mybuffer [MY_IO_BUFSIZE * MY_IO_BUFSIZE];
	out_file.rdbuf()->pubsetbuf(mybuffer,MY_IO_BUFSIZE * MY_IO_BUFSIZE);
	//
#endif
	//========================================================================
	out_file << "GeoSys-OUT: Output ------------------------------------------------\n";
	//========================================================================
	// OUT vector
	size_t out_vector_size (out_vector.size());
	for(size_t i = 0; i < out_vector_size; i++)
		out_vector[i]->Write(&out_file);
	out_file << "#STOP";
	out_file.close();
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   08/2004 OK Implementation
   08/2004 WW Output by steps given in .out file
   03/2005 OK RFO format
   05/2005 OK MSH
   05/2005 OK Profiles at surfaces
   12/2005 OK VAR,MSH,PCS concept
   03/2006 WW Flag to remove existing files
   08/2006 OK FLX calculations
   08/2007 WW Output initial values of variables
**************************************************************************/
void OUTData(double time_current, int time_step_number)
{
	COutput* m_out = NULL;
	CRFProcess* m_pcs = NULL;
	CFEMesh* m_msh = NULL;
	bool OutputBySteps = false;
	double tim_value;

	for (size_t i = 0; i < out_vector.size(); i++)
	{
		m_out = out_vector[i];
		// MSH
		//		m_msh = m_out->GetMSH();
		m_msh = m_out->getMesh();
		if (!m_msh)
			cout << "Warning in OUTData - no MSH data" << endl;
		//OK continue;
		// PCS
		if (m_out->_nod_value_vector.size() > 0)
			m_pcs = m_out->GetPCS(m_out->_nod_value_vector[0]);
		if (m_out->getElementValueVector().size() > 0)
			m_pcs = m_out->GetPCS_ELE(m_out->getElementValueVector()[0]);
		if (!m_pcs)
			m_pcs = m_out->GetPCS();  //OK
		if (!m_pcs)
			cout << "Warning in OUTData - no PCS data" << endl;
		//OK4704 continue;
		//--------------------------------------------------------------------
		m_out->setTime (time_current);
		size_t no_times (m_out->time_vector.size());
		//--------------------------------------------------------------------
		if (no_times == 0 && (m_out->nSteps > 0) && (time_step_number
		                                             % m_out->nSteps == 0))
			OutputBySteps = true;
		if (time_step_number == 0) //WW
			OutputBySteps = true;
		//======================================================================
		// TECPLOT
		if (m_out->dat_type_name.compare("TECPLOT") == 0
		    || m_out->dat_type_name.compare("MATLAB") == 0)
		{
			//			m_out->matlab_delim = " ";
			//			if (m_out->dat_type_name.compare("MATLAB") == 0) // JTARON, just for commenting header for matlab
			//				m_out->matlab_delim = "%";

			switch (m_out->getGeoType())
			{
			case GEOLIB::GEODOMAIN: // domain data
				cout << "Data output: Domain" << endl;
				if (OutputBySteps)
				{
					if (m_out->_pcon_value_vector.size() > 0)
						m_out->PCONWriteDOMDataTEC();  //MX
					else
					{
						m_out->NODWriteDOMDataTEC();
						m_out->ELEWriteDOMDataTEC();
					}
					OutputBySteps = false;
					if (!m_out->_new_file_opened)
						//WW
						m_out->_new_file_opened = true;
				}
				else
				{
					for (size_t j = 0; j < no_times; j++)
						if ((time_current > m_out->time_vector[j]) || fabs(
						            time_current - m_out->time_vector[j])
						    < MKleinsteZahl) //WW MKleinsteZahl
						{
							if (m_out->_pcon_value_vector.size() > 0)
								//MX
								m_out->PCONWriteDOMDataTEC();
							else
							{
								m_out->NODWriteDOMDataTEC();
								m_out->ELEWriteDOMDataTEC();
							}
							m_out->time_vector.erase(
							        m_out->time_vector.begin()
							        + j);
							if (!m_out->_new_file_opened)
								//WW
								m_out->_new_file_opened = true;
							break;
						}
				}
				break;
			//------------------------------------------------------------------
			case GEOLIB::POLYLINE: // profiles along polylines
				std::cout << "Data output: Polyline profile - "
				          << m_out->getGeoName() << std::endl;
				if (OutputBySteps)
				{
					tim_value = m_out->NODWritePLYDataTEC(time_step_number);
					if (tim_value > 0.0)
						//OK
						m_out->TIMValue_TEC(tim_value);
					if (!m_out->_new_file_opened)
						//WW
						m_out->_new_file_opened = true;
					OutputBySteps = false;
				}
				else
				{
					for (size_t j = 0; j < no_times; j++)
						if ((time_current > m_out->time_vector[j]) || fabs(
						            time_current - m_out->time_vector[j])
						    < MKleinsteZahl) //WW MKleinsteZahl
						{
							//OK
							tim_value = m_out->NODWritePLYDataTEC(j + 1);
							if (tim_value > 0.0)
								m_out->TIMValue_TEC(tim_value);
							m_out->time_vector.erase(
							        m_out->time_vector.begin()
							        + j);
							if (!m_out->_new_file_opened)
								//WW
								m_out->_new_file_opened = true;
							break;
						}
				}
				//..............................................................
				break;
			//------------------------------------------------------------------
			case GEOLIB::POINT: // breakthrough curves in points
				cout << "Data output: Breakthrough curves - "
				     << m_out->getGeoName() << endl;
				m_out->NODWritePNTDataTEC(time_current, time_step_number);
				if (!m_out->_new_file_opened)
					m_out->_new_file_opened = true;  //WW
				break;
			//------------------------------------------------------------------
			case GEOLIB::SURFACE: // profiles at surfaces
				cout << "Data output: Surface profile" << endl;
				//..............................................................
				//				if (m_out->_dis_type_name.compare("AVERAGE") == 0) {
				if (m_out->getProcessDistributionType() == FiniteElement::AVERAGE)
				{
					if (OutputBySteps)
					{
						m_out->NODWriteSFCAverageDataTEC(time_current,
						                                 time_step_number);
						OutputBySteps = false;
						if (!m_out->_new_file_opened)
							//WW
							m_out->_new_file_opened = true;
					}
				}
				//..............................................................
				else
				{
					if (OutputBySteps)
					{
						m_out->NODWriteSFCDataTEC(time_step_number);
						OutputBySteps = false;
						if (!m_out->_new_file_opened)
							//WW
							m_out->_new_file_opened = true;
					}
					else
					{
						for (size_t j = 0; j < no_times; j++)
							if ((time_current >
							     m_out->time_vector[j]) || fabs(
							            time_current -
							            m_out->time_vector[j])
							    < MKleinsteZahl) //WW MKleinsteZahl                m_out->NODWriteSFCDataTEC(j);
							{
								m_out->NODWriteSFCDataTEC(j);
								m_out->time_vector.erase(
								        m_out->time_vector.begin()
								        + j);
								if (!m_out->_new_file_opened)
									//WW
									m_out->_new_file_opened =
									        true;
								break;
							}
					}
				}
				//..............................................................
				// ELE data
				if (m_out->getElementValueVector().size() > 0)
					m_out->ELEWriteSFC_TEC();
				//..............................................................
				break;

			//			case 'Y': // Layer
			//				cout << "Data output: Layer" << endl;
			//				if (OutputBySteps) {
			//					m_out->NODWriteLAYDataTEC(time_step_number);
			//					OutputBySteps = false;
			//				} else {
			//					for (j = 0; j < no_times; j++) {
			//						if ((time_current > m_out->time_vector[j]) || fabs(
			//								time_current - m_out->time_vector[j])
			//								<MKleinsteZahl) {
			//							m_out->NODWriteLAYDataTEC(j);
			//							m_out->time_vector.erase(m_out->time_vector.begin()
			//									+ j);
			//							break;
			//						}
			//					}
			//				}
			//
			//				break;
			//------------------------------------------------------------------

			default:
				break;
			}
		}
		//--------------------------------------------------------------------
		// vtk
		else if (m_out->dat_type_name.compare("VTK") == 0)
		{
			switch (m_out->getGeoType())
			{
			case GEOLIB::GEODOMAIN: // domain data
				if (OutputBySteps)
				{
					OutputBySteps = false;
					//OK
					//m_out->WriteDataVTK(time_step_number);
					LegacyVtkInterface vtkOutput(m_msh,
					                             m_out->_nod_value_vector,
					                             m_out->_ele_value_vector,
					                             m_out->mmp_value_vector,
					                             m_out->msh_type_name,
					                             m_out);
					vtkOutput.WriteDataVTK(time_step_number,
					                       m_out->_time,
					                       m_out->file_base_name);

					if (!m_out->_new_file_opened)
						//WW
						m_out->_new_file_opened = true;
				}
				else
				{
					for (size_t j = 0; j < no_times; j++)
						if (time_current >= m_out->time_vector[j])
						{
							//OK
							//m_out->WriteDataVTK(time_step_number);
							LegacyVtkInterface vtkOutput(
							        m_msh,
							        m_out->
							        _nod_value_vector,
							        m_out->
							        _ele_value_vector,
							        m_out->
							        mmp_value_vector,
							        m_out->
							        msh_type_name,
							        m_out);
							vtkOutput.WriteDataVTK(
							        time_step_number,
							        m_out->_time,
							        m_out->
							        file_base_name);
							m_out->time_vector.erase(
							        m_out->time_vector.begin()
							        + j);
							if (!m_out->_new_file_opened)
								//WW
								m_out->_new_file_opened = true;
							break;
						}
				}
				break;
			default:
				break;
			}
		}                           // PVD (ParaView)
		else if (m_out->dat_type_name.find("PVD") != string::npos)
		{
			if (m_out->vtk == NULL)
				m_out->vtk = new CVTK();
			CVTK* vtk = m_out->vtk;

			bool vtk_appended = false;
			if (m_out->dat_type_name.find("PVD_A") != string::npos)
				vtk_appended = true;

			stringstream stm;
			string pvd_vtk_file_name; //, vtk_file_name;

			switch (m_out->getGeoType())
			{
			case GEOLIB::GEODOMAIN: // domain data
				//          static CVTK vtk;
				if (time_step_number == 0)
				{
					std::string pcs_type ("");
					if (m_out->getProcessType() != INVALID_PROCESS)
						pcs_type = convertProcessTypeToString (
						        m_out->getProcessType());
					vtk->InitializePVD(m_out->file_base_name,
					                   pcs_type,
					                   vtk_appended);
				}
				// Set VTU file name and path
				//vtk_file_name = m_out->file_base_name;
				//				if (m_out->pcs_type_name.size() > 0) // PCS
				////				if (m_out->getProcessType() != INVALID_PROCESS)
				//					vtk_file_name += "_" + convertProcessTypeToString (m_out->getProcessType());
				pvd_vtk_file_name = vtk->pvd_vtk_file_name_base;
				stm << time_step_number;
				//vtk_file_name += stm.str() + ".vtu";
				pvd_vtk_file_name += stm.str() + ".vtu";
				//vtk_file_name = pvd_vtk_file_name;
				// Output
				if (OutputBySteps)
				{
					OutputBySteps = false;
					vtk->WriteXMLUnstructuredGrid(pvd_vtk_file_name, m_out,
					                              time_step_number);
					VTK_Info dat;
					dat.timestep = m_out->getTime();
					dat.vtk_file = pvd_vtk_file_name;
					vtk->vec_dataset.push_back(dat);
					vtk->UpdatePVD(vtk->pvd_file_name, vtk->vec_dataset);
				}
				else
				{
					for (size_t j = 0; j < no_times; j++)
						if (time_current >= m_out->time_vector[j])
						{
							vtk->WriteXMLUnstructuredGrid(
							        pvd_vtk_file_name,
							        m_out,
							        time_step_number);
							m_out->time_vector.erase(
							        m_out->time_vector.begin()
							        + j);
							VTK_Info dat;
							dat.timestep = m_out->getTime();
							dat.vtk_file = pvd_vtk_file_name;
							vtk->vec_dataset.push_back(dat);
							vtk->UpdatePVD(vtk->pvd_file_name,
							               vtk->vec_dataset);
							break;
						}
				}
				break;

			default:
				break;
			}
		}
		//		// ROCKFLOW
		//		else if (m_out->dat_type_name.compare("ROCKFLOW") == 0) {
		//			switch (m_out->getGeoType()) {
		//			case GEOLIB::GEODOMAIN: // domain data
		//				if (OutputBySteps) {
		//					OutputBySteps = false;
		//					m_out->WriteRFO(); //OK
		//					if (!m_out->_new_file_opened)
		//						m_out->_new_file_opened = true; //WW
		//				} else {
		//					for (size_t j = 0; j < no_times; j++) {
		//						if ((time_current > m_out->time_vector[j]) ||
		//								fabs(time_current - m_out->time_vector[j]) <MKleinsteZahl) { //WW MKleinsteZahl
		//							m_out->WriteRFO(); //OK
		//							m_out->time_vector.erase(m_out->time_vector.begin()
		//									+ j);
		//							if (!m_out->_new_file_opened)
		//								m_out->_new_file_opened = true; //WW
		//							break;
		//						}
		//					}
		//				}
		//				break;
		//
		//			default:
		//				break;
		//			}
		//		}
		//--------------------------------------------------------------------
		// ELE values
		m_out->CalcELEFluxes();
	}                                     // OUT loop
	//======================================================================
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   01/2005 OK Implementation
   last modified:
**************************************************************************/
void OUTDelete()
{
	const size_t no_out = out_vector.size();
	for(size_t i = 0; i < no_out; i++)
		delete out_vector[i];
	out_vector.clear();
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   05/2005 OK Implementation
   last modification: 03/2010 JTARON
   09/2010 TF
**************************************************************************/
COutput* OUTGet(const std::string & out_name)
{
	ProcessType pcs_type (convertProcessType (out_name));
	for (size_t i = 0; i < out_vector.size(); i++)
		if (out_vector[i]->getProcessType() == pcs_type)
			return out_vector[i];
	return NULL;
}

/**************************************************************************
   FEMLib-Method:
   Task: Return output object for variable name (based on OUTGet)
   Programing:
   03/2010 JTARON Implementation
   last modification:
**************************************************************************/
COutput* OUTGetRWPT(const std::string & out_name)
{
	for (size_t i = 0; i < out_vector.size(); i++)
	{
		COutput* out(out_vector[i]);
		for (size_t j = 0; j < out->getRandomWalkParticleTracingValueVector().size(); j++)
			if (out->getRandomWalkParticleTracingValueVector()[j].compare(out_name) ==
			    0)
				return out;
	}
	return NULL;
}

/***************************************************************************************
   Function checking the consistency of the output data as specified in the input file
   This means up to now, that data for missing processes is not written
   05/09	SB
 *****************************************************************************************/
void OUTCheck()
{
	std::cout << "Checking output data " << std::endl;
	// Go through all out objects (#OUTPUT-section in input file)
	for (size_t i = 0; i < out_vector.size(); i++)
		out_vector[i]->checkConsistency();
}

