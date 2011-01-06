#ifndef VTK_INC
#define VTK_INC

#include <string>
#include <vector>
//#include "rf_out_new.h"
#include "MSHEnums.h"

class COutput;

namespace Mesh_Group
{
   class CFEMesh;
}


using namespace std;
using namespace Mesh_Group;

typedef struct
{
   double timestep;
   string vtk_file;
} VTK_Info;

class CVTK
{
   public:
      vector<VTK_Info> vec_dataset;
      string pvd_file_name;
      string pvd_vtk_file_name_base;
      double useBinary;

      enum VTK_XML_DATA_TYPE { Int8, UInt8, Int16, UInt16, Int32, UInt32, Int64, UInt64, Float32, Float64 };

   protected:
      //for binary output
      bool isInitialized;
      VTK_XML_DATA_TYPE type_UChar;
      VTK_XML_DATA_TYPE type_Int;
      VTK_XML_DATA_TYPE type_UInt;
      VTK_XML_DATA_TYPE type_Long;
      VTK_XML_DATA_TYPE type_Double;
      int SIZE_OF_BLOCK_LENGTH_TAG;
      bool isLittleEndian;                        //Endian(byte order)

   public:
      CVTK(void){isInitialized=false;};
      virtual ~CVTK(void){};

   protected:
      //PVD
      bool WriteHeaderOfPVD(fstream &fin);
      bool WriteEndOfPVD(fstream &fin);
      bool WriteDatasetOfPVD(fstream &fin, double timestep, const string &vtkfile);

      //VTU
      void InitializeVTU();
      unsigned char GetVTKCellType(const MshElemType::type ele_type);
      bool WriteDataArrayHeader(fstream &fin, VTK_XML_DATA_TYPE data_type, const string &str_name, int nr_components, const string &str_format, long offset=-1);
      bool WriteDataArrayFooter(fstream &fin);
      inline bool WriteMeshNodes(fstream &fin, bool output_data, CFEMesh *m_msh, long &offset);
      inline bool WriteMeshElementConnectivity(fstream &fin, bool output_data, CFEMesh *m_msh, long &offset, long &sum_ele_components);
      inline bool WriteMeshElementOffset(fstream &fin, bool output_data, CFEMesh *m_msh, long &offset);
      inline bool WriteMeshElementType(fstream &fin, bool output_data, CFEMesh *m_msh, long &offset);
      inline bool WriteNodalValue(fstream &fin, bool output_data, COutput *out, CFEMesh *m_msh, long &offset);
      inline bool WriteElementValue(fstream &fin, bool output_data, COutput *out, CFEMesh *m_msh, long &offset);

      //util
      template <typename T> void write_value_binary(fstream &fin, T val);
      bool IsLittleEndian();

   public:
      //PVD
      bool InitializePVD(const string &file_base_name, const string &pcs_type_name, bool binary=false);
      bool UpdatePVD(const string &pvdfile, const vector<VTK_Info> &vec_vtk);
      bool CreateDirOfPVD(const string &pvdfile);

      //VTU
      bool WriteXMLUnstructuredGrid(const string &vtkfile, COutput *out, const int time_step_number);

};
#endif
