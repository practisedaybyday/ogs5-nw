/**************************************************************************
FEMLib-Object: MAT-MP
Task: MediumProperties
Programing:
01/2004 OK Implementation
**************************************************************************/

#ifndef rfmat_mp_new_INC
#define rfmat_mp_new_INC
  /* Schutz gegen mehrfaches Einfuegen */

// C++ STL
#include <list>
#include <string>
#include <vector>
#include <fstream>
using namespace std;
// PCSLib
#include "rf_pcs.h"

namespace FiniteElement {class CFiniteElementStd;}
using FiniteElement::CFiniteElementStd;
class CMediumProperties
{
 public:
    CFiniteElementStd *Fem_Ele_Std;
 private:
    friend class FiniteElement::CFiniteElementStd; //WW
    // Data base
	void SetMediumPropertiesDefaultsClay(void);//CMCD 9/2004 GeoSys 4
	void SetMediumPropertiesDefaultsSilt(void);//CMCD 9/2004 GeoSys 4
	void SetMediumPropertiesDefaultsSand(void);//CMCD 9/2004 GeoSys 4
	void SetMediumPropertiesDefaultsGravel(void);//CMCD 9/2004 GeoSys 4
	void SetMediumPropertiesDefaultsCrystalline(void);//CMCD 9/2004 GeoSys 4
	void SetMediumPropertiesDefaultsBordenAquifer(void);//CMCD 9/2004 GeoSys 4
    // Porosity
    double PorosityEffectiveStress(long,double );//CMCD 9/2004 GeoSys 4
    double PorosityVolumetricFreeSwellingConstantIonicstrength(long,double,double);
    double PorosityEffectiveConstrainedSwelling(long,double,double, double*); //MX 1/2005
	double PorosityVolumetricFreeSwelling(long,double,double); //MX 1/2005
	double PorosityEffectiveConstrainedSwellingConstantIonicStrength(long,double,double, double*); //MX 1/2005
    // Permeability
    // Permeabilty stress corrector WW
    int permeability_stress_mode;
    double *c_coefficient;
    //
	double PermeabilityPressureFunctionMethod1(long ,double );//CMCD 9/2004 GeoSys 4
	double PermeabilityPressureFunctionMethod2(long ,double );//CMCD 9/2004 GeoSys 4
	double PermeabilityPressureFunctionMethod3(long ,double );//CMCD 9/2004 GeoSys 4
	double PermeabilityPressureFunctionMethod4(long ,double, double );//CMCD 9/2004 GeoSys 4
    friend class CMediumPropertiesGroup;
 public:
  //-------------------------------------------
  // Methods
  CMediumProperties(void);  // constructor
  ~CMediumProperties(void); // destructor
  CMediumProperties* Get(string);
  CMediumProperties* GetDB(string);
  CMediumProperties* GetByGroupNumber(int);
  void Set(string,string,double);
  void SetDB(string,string,double);
  int GetPropertyType(string);
  ios::pos_type Read(ifstream*);
  void Write(fstream*);
  void WriteTecplot(string);
  double* PermeabilityTensor(long index);
#ifdef RFW_FRACTURE
    double RelativePermeability (long index); //RW/CMCD03/06
#endif
  double Porosity(FiniteElement::CElement* assem=NULL);//CMCD 9/2004 GeoSys 4
  double TortuosityFunction(long number,double*gp,double theta, CFiniteElementStd* assem=NULL);//CMCD 9/2004 GeoSys 4
  double NonlinearFlowFunction(long number, double *gp, double theta);//CMCD 9/2004 GeoSys 4
  double PermeabilityPressureFunction(long index,double *gp,double theta);//CMCD 9/2004 GeoSys 4
  double PermeabilitySaturationFunction(long number,double*gp,double theta,int phase);//CMCD 9/2004 GeoSys 4
  double PermeabilityPorosityFunction(long index,double *gp,double theta);//CMCD 9/2004 GeoSys 4
  double KozenyCarman(double k0/*old permeability*/,
                      double n0/*old porosity*/,
					  double n/*new porosity*/);//HS: 11.2008
  double KozenyCarman_normalized(double k0/*old permeability*/,
                                 double n0/*old porosity*/,
								 double n/*new porosity*/);//HS: 11.2008
  void CalStressPermeabilityFactor(double *kfac, const double T = 273.0); //WW
  void CalStressPermeabilityFactor2(double *kfac, const double T = 273.0); //WW
  void CalStressPermeabilityFactor3(double *kfac); //WW
  void CalStressPermeabilityFactor3_Coef(); //WW
  void CalStressPermeabilityFactor4(double *kfac, double);
#ifdef RFW_FRACTURE
    double PermeabilityFracAperture(long index); //RFW 07/2005
    double CalculateFracAperture(CElem* elem, double search_step); //RFW 04/2005
#endif
  double CapillaryPressureFunction(long number,double*gp,double theta,int phase,double saturation);//CMCD 9/2004 GeoSys 4
  double StorageFunction(long number,double*gp,double theta);//CMCD 9/2004 GeoSys 4
  double HeatCapacity(long number,double theta, CFiniteElementStd* assem=NULL);
  double* HeatConductivityTensor(int number); //MX
  double* HeatDispersionTensorNew(int ip);//CMCD
  double* MassDispersionTensorNew(int ip); //CMCD
  double Density(long number,double*gp,double theta); //OK
  double SaturationCapillaryPressureFunction(long number,double*gp,double theta,int phase); //OK4104
  double SaturationPressureDependency(long number,double*gp,double theta); //OK4104
  double SaturationCapillaryPressureFunction(const double capillary_pressure, const int phase); //WW
  double SaturationPressureDependency(double saturation, double density_fluid, double theta); // WW
	double PressureSaturationDependency(double saturation, double density_fluid, double theta);	// PCH
  double PermeabilitySaturationFunction(const double Saturation, int phase); //WW
  double PorosityVolumetricChemicalReaction(long);  //MX 1/2005
  double Porosity(long number,double theta);//CMCD 9/2004 GeoSys 4
  void SetConstantELEarea(double area, int group);//CMCD 4/2005
  void SetDistributedELEProperties(string); //OK
  void WriteTecplotDistributedProperties(); //OK
  //-------------------------------------------
  // Properties
  // PCS
  string pcs_type_name;   //YD
  CRFProcess*m_pcs; //OK
  CRFProcess*m_pcs_tmp; //HS
  vector<string>pcs_name_vector;
  vector<string>porosity_pcs_name_vector;
  //....................................................................
  //GEO
  string geo_type_name;
  string geo_name;
  vector<string>geo_name_vector; //OK
  int geo_dimension;
  double geo_area;
  string geo_area_file; //OK
  CFEMesh* m_msh; //OK
  //....................................................................
  double density;
  string name;
  int number;
  int porosity_model; // porosity
  int porosity_curve;
  double porosity_model_values[15];
  double porosity;
  double KC_porosity_initial; // HS 11.2008
  double KC_permeability_initial; // HS 11.2008
  string porosity_file; //OK/MB
  int tortuosity_model;
  double tortuosity_model_values[10];
  double tortuosity;
  int flowlinearity_model;
  double flowlinearity_model_values[10];
  int storage_model; // storativity
  double storage_model_values[10];
  double storage;
  int conductivity_model;
  double conductivity;
  int unconfined_flow_group;
  int permeability_model; // permeability
  double permeability;
  double permeability_tensor[10]; // permeability_tensor[9] is the minmum permeability. WW
  string permeability_tensor_type_name;
  string tortuosity_tensor_type_name;
  int permeability_tensor_type;
  int tortuosity_tensor_type;

#ifdef RFW_FRACTURE
  vector<string > relative_permeability_function;//RW/CMCD 03/06
#endif

#ifdef RFW_FRACTURE
    //---------------------------  RFW 07/2005
    string frac_perm_average_type;
    string roughness;
   //---------------------------  RFW 11/2005
    long frac_num, fracs_set;
    vector<string> frac_names;
    vector<double> frac_perm;
    vector<double> avg_aperture;
    vector<double> closed_fraction;
    //---------------------------  RFW 11/2005
#endif
  int permeability_pressure_model;
  double permeability_pressure_model_values[10];
  double permeability_pressure_rel;
  int permeability_saturation_model[3];
  double permeability_saturation;
  string permeability_file; //SB //OK/MB string permeability_dis_type_file;
  string tortuosity_file;	// PCH
  int capillary_pressure_model;
  double capillary_pressure;
  int permeability_porosity_model;
  double permeability_porosity_model_values[10];
  double storativity;
  double saturation_res[3]; // saturation max 3 phases
  double saturation_max[3];
  double saturation_eff[3];
  double saturation_exp[3];
  double saturation_alpha[3]; // JOD
  double heat_capacity; // thermal properties
	int mass_dispersion_model;
  double mass_dispersion_longitudinal;
  double mass_dispersion_transverse;
	double lgpn; //local grid peclet number
  int heat_dispersion_model;
  double heat_dispersion_longitudinal;
  double heat_dispersion_transverse;
  double heat_conductivity_tensor[9];
  int fct_number; // functions
  double permeability_saturation_model_values[2];
  double capillary_pressure_model_values[1];
  double permeability_exp[1]; // JOD
  double permeability_alpha[1]; // JOD
  int heat_diffusion_model;
  //aux
  int m_color[3];
  bool selected;
  int mode;
  string het_file_name; //SB
  // surface water
  double friction_coefficient, friction_exp_slope, friction_exp_depth;  // JOD
  double overland_width, rill_height, rill_epsilon; // JOD
  bool channel;
  double argument; //OK
  //Dual Richards transfer coefficient  YD
  double transfer_coefficient;
  //double unsaturated_hydraulic_conductivity;
  double specific_storage;
  int vol_mat_model; // CB
  double vol_mat; //SB
  int vol_bio_model; // CB
  double vol_bio; //SB
  double foc; // organic carbon content
};

 class CMediumPropertiesGroup  //YD
{
  public:
    CMediumPropertiesGroup() {OrigSize=0;}
    void Set(CRFProcess* m_pcs);
    string pcs_name;
    string pcs_type_name;
    CFEMesh* m_msh;
    vector<CMediumProperties*>mmp_group_vector;
  private:
    int OrigSize;  // For excavation simulation.
};

extern CMediumPropertiesGroup* MMPGetGroup(const string &pcs_type_name); //YD
extern list<CMediumPropertiesGroup*>mmp_group_list;   //YD
extern void MMPGroupDelete(/*string pcs_type_name*/);

extern vector<CMediumProperties*>mmp_vector;
extern void MATLoadDB(string);
extern list<string> keywd_list; //keyword-referenzliste "kw"
extern void comp_keywd_list(string);
extern void read_keywd_list(void);
extern list<string>mat_name_list;

extern void MMPWrite(string);
extern bool MMPRead(string);
extern void MMPWriteTecplot(string);
extern void MMPDelete();
extern CMediumProperties* MMPGet(const string&);
extern void MMP2PCSRelation(CRFProcess*);
extern void GetHeterogeneousFields(); //SB
extern long GetNearestHetVal2(long EleIndex, CFEMesh *m_msh, vector <double> xvals,  vector <double> yvals,  vector <double> zvals,  vector <double> mmpvals);
double GetAverageHetVal2(long EleIndex, CFEMesh *m_msh, vector <double> xvals,  vector <double> yvals,  vector <double> zvals,  vector <double> mmpvals);
extern bool MMPExist(ifstream *mmp_file); //OK
extern bool MMPExist(); //OK

#define MMP_FILE_EXTENSION ".mmp"

#endif
