/**
 * \file FEMEnums.h
 * 31/08/2010 KR inital implementation
 *
 */

#ifndef FEMENUMS_H
#define FEMENUMS_H

#include <string>
#include <limits>

/** Types of physical processes supported by OpenGeoSys
 * if you change this enum, make sure you apply the changes to
 * the functions convertPorcessType(), convertProcessTypeToString(),
 *  isFlowProcess() and isDeformationProcess()
 */
enum ProcessType
{
	INVALID_PROCESS = 0,//!< INVALID_PROCESS
	AIR_FLOW,           //!< AIR_FLOW
	/// M process, single/multi-phase flow
	DEFORMATION,        //!< DEFORMATION
	DEFORMATION_DYNAMIC, //!< ...
	/// C process, single/multi-phase flow
	DEFORMATION_FLOW,   //!< DEFORMATION_FLOW
	/// H process, incompressible flow
	GROUNDWATER_FLOW,   //!< GROUNDWATER_FLOW
	/// T process, single/multi-phase flow
	HEAT_TRANSPORT,     //!< HEAT_TRANSPORT
	FLUID_FLOW,
	FLUID_MOMENTUM, // BC only
	FLUX,
	/// H process, incompressible flow
	LIQUID_FLOW,        //!< LIQUID_FLOW
	MASS_TRANSPORT,     //!< MASS_TRANSPORT
	MULTI_PHASE_FLOW,   //!< MULTI_PHASE_FLOW
	NO_PCS,              //!< NO_PCS
	/// H process, incompressible flow
	OVERLAND_FLOW,      //!< OVERLAND_FLOW
	PS_GLOBAL,          //!< PS_GLOBAL
	RANDOM_WALK,        //!< RANDOM_WALK
	/// H process, incompressible flow
	RICHARDS_FLOW,      //!< RICHARDS_FLOW
	/// H2 process, compressible flow
	TWO_PHASE_FLOW     //!< TWO_PHASE_FLOW
};

/**
 * convert the given string into the appropriate enum value
 * @param pcs_type_string string describing a process type
 * @return enum value describing process type
 */
ProcessType convertProcessType ( const std::string& pcs_type_string );

/**
 * convert the given enum value into the appropriate string
 * @param pcs_type process type described by the enum ProcessType
 * @return string describing the process type
 */
std::string convertProcessTypeToString ( ProcessType pcs_type );

/**
 * checks if the given pcs_type variable corresponds to a flow type of the enum ProcessType
 * @param pcs_type value of enum ProcessType
 * @return true if pcs_type describes a flow process, else false
 */
bool isFlowProcess (ProcessType pcs_type);

/**
 * checks if the given pcs_type variable corresponds to a deformation type of the enum ProcessType
 * @param pcs_type value of enum ProcessType
 * @return true if pcs_type describes a deformation process, else false
 */
bool isDeformationProcess (ProcessType pcs_type);

/**
 * \enum contains all values for primary variables actually handled by OGS
 */
enum PrimaryVariable
{
	INVALID_PV	= 0,                                                            //!< INVALID_PV
	/// Flow (phase)
	PRESSURE,                                                                  //!< PRESSURE
	PRESSURE2,                                                                 //!< PRESSURE2
	PRESSURE_RATE1,	// OUT
	SATURATION,                                                                //!< SATURATION
	SATURATION2,                                                               //!< SATURATION2
	/// Heat transport
	TEMPERATURE,                                                               //!< TEMPERATURE
	/// Deformation
	DISPLACEMENT_X,                                                            //!< DISPLACEMENT_X
	/// Deformation
	DISPLACEMENT_Y,                                                            //!< DISPLACEMENT_Y
	/// Deformation
	DISPLACEMENT_Z,                                                            //!< DISPLACEMENT_Z
	/// Mass transport
	CONCENTRATION,                                                             //!< CONCENTRATION
	HEAD,                                                                      //!< HEAD
	VELOCITY_DM_X,                                                             //!< VELOCITY_DM_X
	VELOCITY_DM_Y,                                                             //!< VELOCITY_DM_Y
	VELOCITY_DM_Z,                                                             //!< VELOCITY_DM_Z
	VELOCITY1_X,
	VELOCITY1_Y,
	VELOCITY1_Z,
	STRESS_XX,			// IC
	STRESS_XY,			// IC
	STRESS_XZ,			// IC
	STRESS_YY,			// IC
	STRESS_YZ,			// IC
	STRESS_ZZ,			// IC
	STRAIN_XX,			// Output
	STRAIN_XY,			// Output
	STRAIN_XZ,			// Output
	STRAIN_YY,			// Output
	STRAIN_YZ,			// Output
	STRAIN_ZZ,			// Output
	STRAIN_PLS,			// Output
	ACCELERATION_X1,                                                           //!< ACCELERATION_X1
	ACCELERATION_Y1,                                                           //!< ACCELERATION_Y1
	ACCELERATION_Z1,                                                           //!< ACCELERATION_Z1
	EXCAVATION,	// ST
#ifdef BRNS
	C3H5O3,
	O2,
	NO3,
	SO4,
	Hplus,
	OH,
	H2CO2,
	H2CO3,
	HCO3,
	CO3,
	Mn2,
	H2S,
	HS,
	S2minus,
	C3H6O3,
	C3H5O2,
	C3H6O2,
	C2H3O2,
	C2H4O2,
	FeS,
	FeCO3,
	FeOOH,
	C3H5O3bio,
	C3H5O2bio,
	O2bio,
	NO3bio,
	SO4bio,
	MnO2bio,
	FeOOHbio,
	BaCO,
	BaCn, //!< Barium Copernicum or Barium Carbon Nitrogen ?
	BaCC,
	BaCm, //!< Barium Curium
	BaCf, //!< Barium Californium
	BaCs, //!< Barium Caesium
	MnO2,
	CaCO3,
	tracer_mob,
#endif
	Isochlor,		                                                           //!< Isochlor
	Decay, //!< Decay used in boundary conditions (only in chemical benchmarks?)
	SorbLin, //!< SorbLin used in boundary conditions (only in chemical benchmarks?)
	SorbDecay,
	SorbFr,
	DecaySorbnonlin, //!< DecaySorbnonlin used in boundary conditions (only in chemical benchmarks?)
	ConsTracer, //!< ConsTracer used in boundary conditions (only in chemical benchmarks?)
	PCE_l, //!< PCE_l used in boundary conditions (only in chemical benchmarks?)
	PCE_h, //!< PCE_h used in boundary conditions (only in chemical benchmarks?)
	TCE_l, //!< TPCE_l used in boundary conditions (only in chemical benchmarks?)
	TCE_h, //!< TCE_h used in boundary conditions (only in chemical benchmarks?)
	TCE_mob, //!< TCE_mob used in boundary conditions (only in chemical benchmarks?)
	TCE_mat, //!< TCE_mob used in boundary conditions (only in chemical benchmarks?)
	DCE_l, //!< DCE_l used in boundary conditions (only in chemical benchmarks?)
	DCE_h, //!< DCE_h used in boundary conditions (only in chemical benchmarks?)
	VC_l, //!< VC_l used in boundary conditions (only in chemical benchmarks?)
	VC_h, //!< VC_h used in boundary conditions (only in chemical benchmarks?)
	Eth_l, //!< Eth_l used in boundary conditions (only in chemical benchmarks?)
	Eth_h, //!< Eth_h used in boundary conditions (only in chemical benchmarks?)
	Micros,
	Xylene,
	Oxygen,
	CO2,
	Fe2,
	Sulfate,
	S,
	Metab_O,
	Biomass,
	Biomass_Aerobic,
	Biomass_Ironreducer,
	Goethite,
	bioFe3,
	Metab_Fe3,
	Meta_ProductS,
	Biomass_Sulfatereducer,
	Toluene,
	Nitrate,
	Ca,
	CaX2,
	Cl,
	Comp_A,
	Comp_B,
	Comp_C,
	C_Cl,
	C_Cs,
	C_HTO,
	C_Na,
	Aromaticum,
	Reaction_Mark,
	TRACER,		// BC
	TRACER1,		// BC
	TRACER_BR,	// BC
	Cons_Tracer,
	trans_DCE_mob,
	trans_DCE_mat,
	cis_DCE_mob,
	cis_DCE_mat,
	Elf_DCE_mob,
	Elf_DCE_mat,
	Acet_mob,
	Acet_mat,
	ClAcet_mob,
	ClAcet_mat,
	VC_mob,
	VC_mat,
	Ethylen_mob,
	Ethylen_mat,
	C4_mob,
	Ethan_mob,
	mobile,
	Mg,
	C_brace_4_brace,
	pH,
	pe,
	Calcite,
	Dolomite_dis,
	K,
	KX,
	N_brace_5_brace,
	Na,
	NaX,
	one_minus_C,
	two_minus_Ca,
	three_minus_Cl,
	four_minus_H,
	five_minus_Mg,
	six_minus_O,
	seven_minus_Sn,
	eight_minus_Zz
};

/**
 * Converts the given string into the appropriate enum value.
 * @param pcs_pv_string string describing the primary variable
 * @return enum value describing the primary variable of the process
 */
PrimaryVariable convertPrimaryVariable ( const std::string& pcs_pv_string );//!< PrimaryVariable

/**
 * Converts the given enum value into the appropriate string.
 * @param pcs_pv primary variable described by the enum ProcessType
 * @return string describing the process type
 */
std::string convertPrimaryVariableToString ( PrimaryVariable pcs_pv );

enum DistributionType
{
	INVALID_DIS_TYPE = 0,
	ANALYTICAL,			// ST
	AVERAGE,
	CONSTANT,			// IC, BC, ST
	CONSTANT_GEO,
	GRADIENT,			// IC
	RESTART,			// IC
	LINEAR,				// BC, ST
	POINT,				// BC
	CONSTANT_NEUMANN,	// ST
	LINEAR_NEUMANN,		// ST
	NORMALDEPTH,		// ST
	CRITICALDEPTH,		// ST
	GREEN_AMPT,			// ST
	SYSTEM_DEPENDENT	// ST
};

/**
 * Converts the given string into the appropriate enum value.
 * @param pcs_pv_string string describing the primary variable
 * @return enum value describing the primary variable of the process
 */
DistributionType convertDisType ( const std::string& dis_type_string );

/**
 * Converts the given enum value into the appropriate string.
 * @param pcs_pv primary variable described by the enum ProcessType
 * @return string describing the process type
 */
std::string convertDisTypeToString ( DistributionType dis_type );


#endif //FEMENUMS_H
