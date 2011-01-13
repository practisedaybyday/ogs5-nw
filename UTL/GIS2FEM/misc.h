#ifndef misc_INC
#define misc_INC
#include <vector>
using namespace std;

namespace  Mesh_Group{class CFEMesh; }
using  Mesh_Group::CFEMesh;
extern vector<Mesh_Group::CFEMesh*> fem_msh_vector;

extern bool FEMRead(string);
#endif