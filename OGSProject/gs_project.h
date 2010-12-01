/**************************************************************************
GeoSys-Project
Task: class implementation
Programing:
01/2005 OK Implementation
last modified:
**************************************************************************/
#ifndef gs_project_INC
#define gs_project_INC

#include "Configure.h"
// C++ STL
#include <string>
#include <vector>
using namespace std;
// GEOLib
// FEMLib
//----------------------------------------------------------------
// ToDo CGSProject
class CGSProject
{
  private:
  public:
    string path;
    string base;
    string type;
};

extern vector<CGSProject*>gsp_vector;
extern void GSPRead(string);
extern void GSPWrite(); 
extern void GSPAddMember(string);
extern void GSPRemoveMember(string);
extern void GSPWriteData();
extern CGSProject* GSPGetMember(string);
extern void GSPAddMemberNew(string path_base_orig,string path_base_copy,string type);
extern bool GSPSimulatorReady();
extern bool GSPReadData(string); //OK

extern string g_gsp_path;
extern string g_gsp_base;
#define GSP_FILE_EXTENSION ".gsp"

#endif
