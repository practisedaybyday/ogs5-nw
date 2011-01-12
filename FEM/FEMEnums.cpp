/*
 * FEMEnums.cpp
 *
 *  Created on: Sep 2, 2010
 *      Author: TF
 */

#include "FEMEnums.h"
#include <iostream>
#include <cstdlib>

ProcessType convertProcessType ( const std::string& pcs_type_string )
{
   if (pcs_type_string.compare ("LIQUID_FLOW") == 0) return LIQUID_FLOW;
   if (pcs_type_string.compare ("FLUID_FLOW") == 0) return FLUID_FLOW;
   if (pcs_type_string.compare ("TWO_PHASE_FLOW") == 0) return TWO_PHASE_FLOW;
   if (pcs_type_string.compare ("RICHARDS_FLOW") == 0) return RICHARDS_FLOW;
   if (pcs_type_string.compare ("OVERLAND_FLOW") == 0) return OVERLAND_FLOW;
   if (pcs_type_string.compare ("GROUNDWATER_FLOW") == 0) return GROUNDWATER_FLOW;
   if (pcs_type_string.compare ("HEAT_TRANSPORT") == 0) return HEAT_TRANSPORT;
   if (pcs_type_string.compare ("DEFORMATION") == 0) return DEFORMATION;
   if (pcs_type_string.compare ("DEFORMATION_FLOW") == 0) return DEFORMATION_FLOW;
   if (pcs_type_string.compare ("DEFORMATION_DYNAMIC") == 0) return DEFORMATION_DYNAMIC;
   if (pcs_type_string.compare ("MASS_TRANSPORT") == 0) return MASS_TRANSPORT;
   if (pcs_type_string.compare ("MULTI_PHASE_FLOW") == 0) return MULTI_PHASE_FLOW;
   if (pcs_type_string.compare ("AIR_FLOW") == 0) return AIR_FLOW;
   if (pcs_type_string.compare ("FLUID_MOMENTUM") == 0) return FLUID_MOMENTUM;
   if (pcs_type_string.compare ("RANDOM_WALK") == 0) return RANDOM_WALK;
   if (pcs_type_string.compare ("FLUX") == 0) return FLUX;
   if (pcs_type_string.compare ("PS_GLOBAL") == 0) return PS_GLOBAL;
   if (pcs_type_string.compare ("NO_PCS") == 0) return NO_PCS;
   else
   {
      std::cout << "WARNING in convertProcessType: process type #" << pcs_type_string << "# unknown" << std::endl;
      exit (1);
   }
   return INVALID_PROCESS;
}


std::string convertProcessTypeToString ( ProcessType pcs_type )
{
   if (pcs_type == LIQUID_FLOW) return "LIQUID_FLOW";
   if (pcs_type == FLUID_FLOW) return "FLUID_FLOW";
   if (pcs_type == TWO_PHASE_FLOW) return "TWO_PHASE_FLOW";
   if (pcs_type == RICHARDS_FLOW) return "RICHARDS_FLOW";
   if (pcs_type == OVERLAND_FLOW) return "OVERLAND_FLOW";
   if (pcs_type == GROUNDWATER_FLOW) return "GROUNDWATER_FLOW";
   if (pcs_type == HEAT_TRANSPORT) return "HEAT_TRANSPORT";
   if (pcs_type == DEFORMATION) return "DEFORMATION";
   if (pcs_type == DEFORMATION_FLOW) return "DEFORMATION_FLOW";
   if (pcs_type == DEFORMATION_DYNAMIC) return "DEFORMATION_DYNAMIC";
   if (pcs_type == MASS_TRANSPORT) return "MASS_TRANSPORT";
   if (pcs_type == MULTI_PHASE_FLOW) return "MULTI_PHASE_FLOW";
   if (pcs_type == AIR_FLOW) return "AIR_FLOW";
   if (pcs_type == FLUID_MOMENTUM) return "FLUID_MOMENTUM";
   if (pcs_type == RANDOM_WALK) return "RANDOM_WALK";
   if (pcs_type == FLUX) return "FLUX";
   if (pcs_type ==   PS_GLOBAL) return "PS_GLOBAL";
   if (pcs_type ==   NO_PCS) return "NO_PCS";
   return "INVALID_PROCESS";
}


bool isFlowProcess (ProcessType pcs_type)
{
   if (pcs_type == LIQUID_FLOW || pcs_type == FLUID_FLOW
      || pcs_type == TWO_PHASE_FLOW || pcs_type == RICHARDS_FLOW
      || pcs_type == OVERLAND_FLOW || pcs_type == GROUNDWATER_FLOW
      || pcs_type == DEFORMATION_FLOW || pcs_type == MULTI_PHASE_FLOW
      || pcs_type == AIR_FLOW)
      return true;
   return false;
}


bool isDeformationProcess (ProcessType pcs_type)
{
   if (pcs_type == DEFORMATION || pcs_type == DEFORMATION_FLOW || pcs_type == DEFORMATION_DYNAMIC)
      return true;
   return false;
}


PrimaryVariable convertPrimaryVariable ( const std::string& pcs_pv_string )
{
   if (pcs_pv_string.compare ("PRESSURE1") == 0) return PRESSURE;
   if (pcs_pv_string.compare ("PRESSURE2") == 0) return PRESSURE2;
   if (pcs_pv_string.compare ("PRESSURE_RATE1") == 0) return PRESSURE_RATE1;
   if (pcs_pv_string.compare ("SATURATION1") == 0) return SATURATION;
   if (pcs_pv_string.compare ("SATURATION2") == 0) return SATURATION2;
   if (pcs_pv_string.compare ("TEMPERATURE1") == 0) return TEMPERATURE;
   if (pcs_pv_string.compare ("DISPLACEMENT_X1") == 0) return DISPLACEMENT_X;
   if (pcs_pv_string.compare ("DISPLACEMENT_Y1") == 0) return DISPLACEMENT_Y;
   if (pcs_pv_string.compare ("DISPLACEMENT_Z1") == 0) return DISPLACEMENT_Z;
   if (pcs_pv_string.compare ("CONCENTRATION1") == 0) return CONCENTRATION;
   if (pcs_pv_string.compare ("HEAD") == 0) return HEAD;
   if (pcs_pv_string.compare ("VELOCITY_DM_X") == 0) return VELOCITY_DM_X;
   if (pcs_pv_string.compare ("VELOCITY_DM_Y") == 0) return VELOCITY_DM_Y;
   if (pcs_pv_string.compare ("VELOCITY_DM_Z") == 0) return VELOCITY_DM_Z;
   if (pcs_pv_string.compare ("VELOCITY1_X") == 0) return VELOCITY1_X;
   if (pcs_pv_string.compare ("VELOCITY1_Y") == 0) return VELOCITY1_Y;
   if (pcs_pv_string.compare ("VELOCITY1_Z") == 0) return VELOCITY1_Z;
   if (pcs_pv_string.compare ("STRESS_XX") == 0) return STRESS_XX;
   if (pcs_pv_string.compare ("STRESS_XY") == 0) return STRESS_XY;
   if (pcs_pv_string.compare ("STRESS_XZ") == 0) return STRESS_XZ;
   if (pcs_pv_string.compare ("STRESS_YY") == 0) return STRESS_YY;
   if (pcs_pv_string.compare ("STRESS_YZ") == 0) return STRESS_YZ;
   if (pcs_pv_string.compare ("STRESS_ZZ") == 0) return STRESS_ZZ;
   //	if (pcs_pv_string.compare ("STRAIN_XX") == 0) return STRAIN_XX;
   //	if (pcs_pv_string.compare ("STRAIN_XY") == 0) return STRAIN_XY;
   //	if (pcs_pv_string.compare ("STRAIN_XZ") == 0) return STRAIN_XZ;
   //	if (pcs_pv_string.compare ("STRAIN_YY") == 0) return STRAIN_YY;
   //	if (pcs_pv_string.compare ("STRAIN_YZ") == 0) return STRAIN_YZ;
   //	if (pcs_pv_string.compare ("STRAIN_ZZ") == 0) return STRAIN_ZZ;
   //	if (pcs_pv_string.compare ("STRAIN_PLS") == 0) return STRAIN_PLS;
   if (pcs_pv_string.compare ("ACCELERATION_X1") == 0) return ACCELERATION_X1;
   if (pcs_pv_string.compare ("ACCELERATION_Y1") == 0) return ACCELERATION_Y1;
   if (pcs_pv_string.compare ("ACCELERATION_Z1") == 0) return ACCELERATION_Z1;
   if (pcs_pv_string.compare ("EXCAVATION") == 0) return EXCAVATION;
   if (pcs_pv_string.compare ("Isochlor") == 0) return Isochlor;
   if (pcs_pv_string.compare ("Decay") == 0) return Decay;
   if (pcs_pv_string.compare ("SorbLin") == 0) return SorbLin;
   if (pcs_pv_string.compare ("SorbDecay") == 0) return SorbDecay;
   if (pcs_pv_string.compare ("SorbFr") == 0) return SorbFr;
   if (pcs_pv_string.compare ("DecaySorbnonlin") == 0) return DecaySorbnonlin;
   if (pcs_pv_string.compare ("PCE-l") == 0) return PCE_l;
   if (pcs_pv_string.compare ("PCE-h") == 0) return PCE_h;
   if (pcs_pv_string.compare ("TCE-l") == 0) return TCE_l;
   if (pcs_pv_string.compare ("TCE-h") == 0) return TCE_h;
   if (pcs_pv_string.compare ("TCE-mob") == 0) return TCE_mob;
   if (pcs_pv_string.compare ("TCE-mat") == 0) return TCE_mat;
   if (pcs_pv_string.compare ("DCE-l") == 0) return DCE_l;
   if (pcs_pv_string.compare ("DCE-h") == 0) return DCE_h;
   if (pcs_pv_string.compare ("VC-l") == 0) return VC_l;
   if (pcs_pv_string.compare ("VC-h") == 0) return VC_h;
   if (pcs_pv_string.compare ("Eth-l") == 0) return Eth_l;
   if (pcs_pv_string.compare ("Eth-h") == 0) return Eth_h;
   if (pcs_pv_string.compare ("Micros") == 0) return Micros;
   if (pcs_pv_string.compare ("Xylene") == 0) return Xylene;
   if (pcs_pv_string.compare ("Tracer") == 0) return TRACER;
   if (pcs_pv_string.compare ("Tracer1") == 0) return TRACER1;
   if (pcs_pv_string.compare ("Tracer_Br") == 0) return TRACER_BR;
   if (pcs_pv_string.compare ("ConsTracer") == 0) return ConsTracer;
   if (pcs_pv_string.compare ("Cons_Tracer") == 0) return Cons_Tracer;
   if (pcs_pv_string.compare ("Oxygen") == 0) return Oxygen;
   if (pcs_pv_string.compare ("CO2") == 0) return CO2;
   if (pcs_pv_string.compare ("Fe2") == 0) return Fe2;
   if (pcs_pv_string.compare ("Sulfate") == 0) return Sulfate;
   if (pcs_pv_string.compare ("S") == 0) return S;
   if (pcs_pv_string.compare ("Metab_O") == 0) return Metab_O;
   if (pcs_pv_string.compare ("Biomass") == 0) return Biomass;
   if (pcs_pv_string.compare ("Biomass_Aerobic") == 0) return Biomass_Aerobic;
   if (pcs_pv_string.compare ("Biomass_Ironreducer") == 0) return Biomass_Ironreducer;
   if (pcs_pv_string.compare ("Goethite") == 0) return Goethite;
   if (pcs_pv_string.compare ("bioFe3") == 0) return bioFe3;
   if (pcs_pv_string.compare ("Metab_Fe3") == 0) return Metab_Fe3;
   if (pcs_pv_string.compare ("Meta_ProductS") == 0) return Meta_ProductS;
   if (pcs_pv_string.compare ("Biomass_Sulfatereducer") == 0) return Biomass_Sulfatereducer;
   if (pcs_pv_string.compare ("Toluene") == 0) return Toluene;
   if (pcs_pv_string.compare ("Nitrate") == 0) return Nitrate;
   if (pcs_pv_string.compare ("Ca") == 0) return Ca;
   if (pcs_pv_string.compare ("ca") == 0) return Ca;
   if (pcs_pv_string.compare ("CaX2") == 0) return CaX2;
   if (pcs_pv_string.compare ("Cl") == 0) return Cl;
   if (pcs_pv_string.compare ("Comp_A") == 0) return Comp_A;
   if (pcs_pv_string.compare ("Comp_B") == 0) return Comp_B;
   if (pcs_pv_string.compare ("Comp_C") == 0) return Comp_C;
   if (pcs_pv_string.compare ("C_Cl") == 0) return C_Cl;
   if (pcs_pv_string.compare ("C_Cs") == 0) return C_Cs;
   if (pcs_pv_string.compare ("C_HTO") == 0) return C_HTO;
   if (pcs_pv_string.compare ("C_Na") == 0) return C_Na;
   if (pcs_pv_string.compare ("Aromaticum") == 0) return Aromaticum;
   if (pcs_pv_string.compare ("Reaction_Mark") == 0) return Reaction_Mark;
   if (pcs_pv_string.compare ("trans-DCE-mob") == 0) return trans_DCE_mob;
   if (pcs_pv_string.compare ("trans-DCE-mat") == 0) return trans_DCE_mat;
   if (pcs_pv_string.compare ("cis-DCE-mob") == 0) return cis_DCE_mob;
   if (pcs_pv_string.compare ("cis-DCE-mat") == 0) return cis_DCE_mat;
   if (pcs_pv_string.compare ("11-DCE-mob") == 0) return Elf_DCE_mob;
   if (pcs_pv_string.compare ("11-DCE-mat") == 0) return Elf_DCE_mat;
   if (pcs_pv_string.compare ("Acet-mob") == 0) return Acet_mob;
   if (pcs_pv_string.compare ("Acet-mat") == 0) return Acet_mat;
   if (pcs_pv_string.compare ("ClAcet-mob") == 0) return ClAcet_mob;
   if (pcs_pv_string.compare ("ClAcet-mat") == 0) return ClAcet_mat;
   if (pcs_pv_string.compare ("VC-mob") == 0) return VC_mob;
   if (pcs_pv_string.compare ("VC-mat") == 0) return VC_mat;
   if (pcs_pv_string.compare ("Ethylen-mob") == 0) return Ethylen_mob;
   if (pcs_pv_string.compare ("Ethylen-mat") == 0) return Ethylen_mat;
   if (pcs_pv_string.compare ("C4-mob") == 0) return C4_mob;
   if (pcs_pv_string.compare ("Ethan-mob") == 0) return Ethan_mob;
   if (pcs_pv_string.compare ("mobile") == 0) return mobile;
   if (pcs_pv_string.compare ("C(4)") == 0) return C_brace_4_brace;
   if (pcs_pv_string.compare ("Mg") == 0) return Mg;
   if (pcs_pv_string.compare ("pH") == 0) return pH;
   if (pcs_pv_string.compare ("pe") == 0) return pe;
   if (pcs_pv_string.compare ("Calcite") == 0) return Calcite;
   if (pcs_pv_string.compare ("Dolomite(dis)") == 0) return Dolomite_dis;
   if (pcs_pv_string.compare ("N(5)") == 0) return N_brace_5_brace;
   if (pcs_pv_string.compare ("K") == 0) return K;
   if (pcs_pv_string.compare ("KX") == 0) return KX;
   if (pcs_pv_string.compare ("Na") == 0) return Na;
   if (pcs_pv_string.compare ("NaX") == 0) return NaX;
   if (pcs_pv_string.compare ("1-C") == 0) return one_minus_C;
   if (pcs_pv_string.compare ("2-Ca") == 0) return two_minus_Ca;
   if (pcs_pv_string.compare ("3-Cl") == 0) return three_minus_Cl;
   if (pcs_pv_string.compare ("4-H") == 0) return four_minus_H;
   if (pcs_pv_string.compare ("5-Mg") == 0) return five_minus_Mg;
   if (pcs_pv_string.compare ("6-O") == 0) return six_minus_O;
   if (pcs_pv_string.compare ("7-Sn") == 0) return seven_minus_Sn;
   if (pcs_pv_string.compare ("8-Zz") == 0) return eight_minus_Zz;
#ifdef BRNS
   if (pcs_pv_string.compare ("c3h5o3") == 0) return C3H5O3;
   if (pcs_pv_string.compare ("o2") == 0) return O2;
   if (pcs_pv_string.compare ("no3") == 0) return NO3;
   if (pcs_pv_string.compare ("so4") == 0) return SO4;
   if (pcs_pv_string.compare ("fe2") == 0) return Fe2;
   if (pcs_pv_string.compare ("hplus") == 0) return Hplus;
   if (pcs_pv_string.compare ("oh") == 0) return OH;
   if (pcs_pv_string.compare ("h2co2") == 0) return H2CO2;
   if (pcs_pv_string.compare ("h2co3") == 0) return H2CO3;
   if (pcs_pv_string.compare ("hco3") == 0) return HCO3;
   if (pcs_pv_string.compare ("co3") == 0) return CO3;
   if (pcs_pv_string.compare ("mn2") == 0) return Mn2;
   if (pcs_pv_string.compare ("h2s") == 0) return H2S;
   if (pcs_pv_string.compare ("hs") == 0) return HS;
   if (pcs_pv_string.compare ("s2minus") == 0) return S2minus;
   if (pcs_pv_string.compare ("c3h6o3") == 0) return C3H6O3;
   if (pcs_pv_string.compare ("c3h5o2") == 0) return C3H5O2;
   if (pcs_pv_string.compare ("c3h6o2") == 0) return C3H6O2;
   if (pcs_pv_string.compare ("c2h3o2") == 0) return C2H3O2;
   if (pcs_pv_string.compare ("c2h4o2") == 0) return C2H4O2;
   if (pcs_pv_string.compare ("fes") == 0) return FeS;
   if (pcs_pv_string.compare ("feco3") == 0) return FeCO3;
   if (pcs_pv_string.compare ("feooh") == 0) return FeOOH;
   if (pcs_pv_string.compare ("c3h5o3bio") == 0) return C3H5O3bio;
   if (pcs_pv_string.compare ("c3h5o2bio") == 0) return C3H5O2bio;
   if (pcs_pv_string.compare ("o2bio") == 0) return O2bio;
   if (pcs_pv_string.compare ("no3bio") == 0) return NO3bio;
   if (pcs_pv_string.compare ("so4bio") == 0) return SO4bio;
   if (pcs_pv_string.compare ("mno2bio") == 0) return MnO2bio;
   if (pcs_pv_string.compare ("feoohbio") == 0) return FeOOHbio;
   if (pcs_pv_string.compare ("baco") == 0) return BaCO;
                                                  //!< Barium Copernicum or Barium Carbon Nitrogen ?
   if (pcs_pv_string.compare ("bacn") == 0) return BaCn;
   if (pcs_pv_string.compare ("bacc") == 0) return BaCC;
                                                  //!< Barium Curium
   if (pcs_pv_string.compare ("bacm") == 0) return BaCm;
                                                  //!< Barium Californium
   if (pcs_pv_string.compare ("bacf") == 0) return BaCf;
                                                  //!< Barium Caesium
   if (pcs_pv_string.compare ("bacs") == 0) return BaCs;
   if (pcs_pv_string.compare ("mno2") == 0) return MnO2;
   if (pcs_pv_string.compare ("caco3") == 0) return CaCO3;
   if (pcs_pv_string.compare ("tracer_mob") == 0) return tracer_mob;
#endif
   else
   {
      std::cout << "convertPrimaryVariable #" << pcs_pv_string << "# not found" << std::endl;
      exit (1);
   }
   return INVALID_PV;
}


std::string convertPrimaryVariableToString ( PrimaryVariable pcs_pv )
{
   if (pcs_pv == PRESSURE) return "PRESSURE1";
   if (pcs_pv == PRESSURE2) return "PRESSURE2";
   if (pcs_pv == PRESSURE_RATE1) return "PRESSURE_RATE1";
   if (pcs_pv == SATURATION) return "SATURATION1";
   if (pcs_pv == SATURATION2) return "SATURATION2";
   if (pcs_pv == TEMPERATURE) return "TEMPERATURE1";
   if (pcs_pv == DISPLACEMENT_X) return "DISPLACEMENT_X1";
   if (pcs_pv == DISPLACEMENT_Y) return "DISPLACEMENT_Y1";
   if (pcs_pv == DISPLACEMENT_Z) return "DISPLACEMENT_Z1";
   if (pcs_pv == CONCENTRATION) return "CONCENTRATION1";
   if (pcs_pv == HEAD) return "HEAD";
   if (pcs_pv == VELOCITY_DM_X) return "VELOCITY_DM_X";
   if (pcs_pv == VELOCITY_DM_Y) return "VELOCITY_DM_Y";
   if (pcs_pv == VELOCITY_DM_Z) return "VELOCITY_DM_Z";
   if (pcs_pv == VELOCITY1_X) return "VELOCITY1_X";
   if (pcs_pv == VELOCITY1_Y) return "VELOCITY1_Y";
   if (pcs_pv == VELOCITY1_Z) return "VELOCITY1_Z";
   if (pcs_pv == STRESS_XX) return "STRESS_XX";
   if (pcs_pv == STRESS_XY) return "STRESS_XY";
   if (pcs_pv == STRESS_XZ) return "STRESS_XZ";
   if (pcs_pv == STRESS_YY) return "STRESS_YY";
   if (pcs_pv == STRESS_YZ) return "STRESS_YZ";
   if (pcs_pv == STRESS_ZZ) return "STRESS_ZZ";
   //	if (pcs_pv == STRAIN_XX) return "STRAIN_XX";
   //	if (pcs_pv == STRAIN_XY) return "STRAIN_XY";
   //	if (pcs_pv == STRAIN_XZ) return "STRAIN_XZ";
   //	if (pcs_pv == STRAIN_YY) return "STRAIN_YY";
   //	if (pcs_pv == STRAIN_YZ) return "STRAIN_YZ";
   //	if (pcs_pv == STRAIN_ZZ) return "STRAIN_ZZ";
   //	if (pcs_pv == STRAIN_PLS) return "STRAIN_PLS";
   if (pcs_pv == ACCELERATION_X1) return "ACCELERATION_X1";
   if (pcs_pv == ACCELERATION_Y1) return "ACCELERATION_Y1";
   if (pcs_pv == ACCELERATION_Z1) return "ACCELERATION_Z1";
   if (pcs_pv == EXCAVATION) return "EXCAVATION";
#ifdef BRNS
   if (pcs_pv == C3H5O3) return "c3h5o3";
   if (pcs_pv == O2) return "o2";
   if (pcs_pv == NO3) return "no3";
   if (pcs_pv == SO4) return "so4";
   if (pcs_pv == Fe2) return "fe2";
   if (pcs_pv == Hplus) return "hplus";
   if (pcs_pv == OH) return "oh";
   if (pcs_pv == H2CO2) return "h2co2";
   if (pcs_pv == H2CO3) return "h2co3";
   if (pcs_pv == HCO3) return "hco3";
   if (pcs_pv == CO3) return "co3";
   if (pcs_pv == Mn2) return "mn2";
   if (pcs_pv == H2S) return "h2s";
   if (pcs_pv == HS) return "hs";
   if (pcs_pv == S2minus) return "s2minus";
   if (pcs_pv == C3H6O3) return "c3h6o3";
   if (pcs_pv == C3H5O2) return "c3h5o2";
   if (pcs_pv == C3H6O2) return "c3h6o2";
   if (pcs_pv == C2H3O2) return "c2h3o2";
   if (pcs_pv == C2H4O2) return "c2h4o2";
   if (pcs_pv == FeS) return "fes";
   if (pcs_pv == FeCO3) return "feco3";
   if (pcs_pv == FeOOH) return "feooh";
   if (pcs_pv == C3H5O3bio) return "c3h5o3bio";
   if (pcs_pv == C3H5O2bio) return "3h5o2bio";
   if (pcs_pv == O2bio) return "o2bio";
   if (pcs_pv == NO3bio) return "no3bio";
   if (pcs_pv == SO4bio) return "so4bio";
   if (pcs_pv == MnO2bio) return "mno2bio";
   if (pcs_pv == FeOOHbio) return "feoohbio";
   if (pcs_pv == BaCO) return "baco";
   if (pcs_pv == BaCn) return "bacn";
   if (pcs_pv == BaCC) return "bacc";
   if (pcs_pv == BaCm) return "bacm";
   if (pcs_pv == BaCf) return "bacf";
   if (pcs_pv == BaCs) return "bacs";
   if (pcs_pv == MnO2) return "mno2";
   if (pcs_pv == CaCO3) return "caco3";
   if (pcs_pv == tracer_mob) return "tracer_mob";
#endif
   if (pcs_pv == Isochlor) return "Isochlor";
   if (pcs_pv == Decay) return "Decay";
   if (pcs_pv == SorbLin) return "SorbLin";
   if (pcs_pv == SorbDecay) return "SorbDecay";
   if (pcs_pv == SorbFr) return "SorbFr";
   if (pcs_pv == DecaySorbnonlin) return "DecaySorbnonlin";
   if (pcs_pv == PCE_l) return "PCE-l";
   if (pcs_pv == PCE_h) return "PCE-h";
   if (pcs_pv == TCE_l) return "TCE-l";
   if (pcs_pv == TCE_h) return "TCE-h";
   if (pcs_pv == TCE_mob) return "TCE-mob";
   if (pcs_pv == TCE_mat) return "TCE-mat";
   if (pcs_pv == DCE_l) return "DCE-l";
   if (pcs_pv == DCE_h) return "DCE-h";
   if (pcs_pv == VC_l) return "VC-l";
   if (pcs_pv == VC_h) return "VC-h";
   if (pcs_pv == Eth_l) return "Eth-l";
   if (pcs_pv == Eth_h) return "Eth-h";
   if (pcs_pv == Micros) return "Micros";
   if (pcs_pv == Xylene) return "Xylene";
   if (pcs_pv == TRACER) return "Tracer";
   if (pcs_pv == TRACER1) return "Tracer1";
   if (pcs_pv == TRACER_BR) return "Tracer_Br";
   if (pcs_pv == ConsTracer) return "ConsTracer";
   if (pcs_pv == Cons_Tracer) return "Cons_Tracer";
   if (pcs_pv == Oxygen) return "Oxygen";
   if (pcs_pv == CO2) return "CO2";
   if (pcs_pv == Fe2) return "Fe2";
   if (pcs_pv == Sulfate) return "Sulfate";
   if (pcs_pv == S) return "S";
   if (pcs_pv == Metab_O) return "Metab_O";
   if (pcs_pv == Biomass) return "Biomass";
   if (pcs_pv == Biomass_Aerobic) return "Biomass_Aerobic";
   if (pcs_pv == Biomass_Ironreducer) return "Biomass_Ironreducer";
   if (pcs_pv == Goethite) return "Goethite";
   if (pcs_pv == bioFe3) return "bioFe3";
   if (pcs_pv == Metab_Fe3) return "Metab_Fe3";
   if (pcs_pv == Meta_ProductS) return "Meta_ProductS";
   if (pcs_pv == Biomass_Sulfatereducer) return "Biomass_Sulfatereducer";
   if (pcs_pv == Toluene) return "Toluene";
   if (pcs_pv == Nitrate) return "Nitrate";
   if (pcs_pv == Ca) return "Ca";
   if (pcs_pv == CaX2) return "CaX2";
   if (pcs_pv == Cl) return "Cl";
   if (pcs_pv == Comp_A) return "Comp_A";
   if (pcs_pv == Comp_B) return "Comp_B";
   if (pcs_pv == Comp_C) return "Comp_C";
   if (pcs_pv == C_Cl) return "C_Cl";
   if (pcs_pv == C_Cs) return "C_Cs";
   if (pcs_pv == C_HTO) return "C_HTO";
   if (pcs_pv == C_Na) return "C_Na";
   if (pcs_pv == Aromaticum) return "Aromaticum";
   if (pcs_pv == Reaction_Mark) return "Reaction_Mark";
   if (pcs_pv == trans_DCE_mob) return "trans-DCE-mob";
   if (pcs_pv == trans_DCE_mat) return "trans-DCE-mat";
   if (pcs_pv == cis_DCE_mob) return "cis-DCE-mob";
   if (pcs_pv == cis_DCE_mat) return "cis-DCE-mat";
   if (pcs_pv == Elf_DCE_mob) return "11-DCE-mob";
   if (pcs_pv == Elf_DCE_mat) return "11-DCE-mat";
   if (pcs_pv == Acet_mob) return "Acet-mob";
   if (pcs_pv == Acet_mat) return "Acet-mat";
   if (pcs_pv == ClAcet_mob) return "ClAcet-mob";
   if (pcs_pv == ClAcet_mat) return "ClAcet-mat";
   if (pcs_pv == VC_mob) return "VC-mob";
   if (pcs_pv == VC_mat) return "VC-mat";
   if (pcs_pv == Ethylen_mob) return "Ethylen-mob";
   if (pcs_pv == Ethylen_mat) return "Ethylen-mat";
   if (pcs_pv == C4_mob) return "C4-mob";
   if (pcs_pv == Ethan_mob) return "Ethan-mob";
   if (pcs_pv == mobile) return "mobile";
   if (pcs_pv == C_brace_4_brace) return "C(4)";
   if (pcs_pv == Mg) return "Mg";
   if (pcs_pv == pH) return "pH";
   if (pcs_pv == pe) return "pe";
   if (pcs_pv == Calcite) return "Calcite";
   if (pcs_pv == Dolomite_dis) return "Dolomite(dis)";
   if (pcs_pv == N_brace_5_brace) return "N(5)";
   if (pcs_pv == K) return "K";
   if (pcs_pv == KX) return "KX";
   if (pcs_pv == Na) return "Na";
   if (pcs_pv == NaX) return "NaX";
   if (pcs_pv == one_minus_C) return "1-C";
   if (pcs_pv == two_minus_Ca) return "2-Ca";
   if (pcs_pv == three_minus_Cl) return "3-Cl";
   if (pcs_pv == four_minus_H) return "4-H";
   if (pcs_pv == five_minus_Mg) return "5-Mg";
   if (pcs_pv == six_minus_O) return "6-O";
   if (pcs_pv == seven_minus_Sn) return "7-Sn";
   if (pcs_pv == eight_minus_Zz) return "8-Zz";
   return "INVALID_PRIMARY_VARIABLE";
}

namespace FiniteElement {

DistributionType convertDisType(const std::string& dis_type_string)
{
	if (dis_type_string.compare("CONSTANT") == 0) return CONSTANT;
	if (dis_type_string.compare("ANALYTICAL") == 0) return ANALYTICAL;
	if (dis_type_string.compare("AVERAGE") == 0) return AVERAGE;
	if (dis_type_string.compare("CONSTANT_GEO") == 0) return CONSTANT_GEO;
	if (dis_type_string.compare("GRADIENT") == 0) return GRADIENT;
	if (dis_type_string.compare("RESTART") == 0) return RESTART;
	if (dis_type_string.compare("LINEAR") == 0) return LINEAR;
	if (dis_type_string.compare("POINT") == 0) return POINT;
	if (dis_type_string.compare("CONSTANT_NEUMANN") == 0) return CONSTANT_NEUMANN;
	if (dis_type_string.compare("LINEAR_NEUMANN") == 0) return LINEAR_NEUMANN;
	if (dis_type_string.compare("NORMALDEPTH") == 0) return NORMALDEPTH;
	if (dis_type_string.compare("CRITICALDEPTH") == 0) return CRITICALDEPTH;
	if (dis_type_string.compare("GREEN_AMPT") == 0) return GREEN_AMPT;
	if (dis_type_string.compare("SYSTEM_DEPENDENT") == 0) return SYSTEM_DEPENDENT;
	if (dis_type_string.compare("PRECIPITATION") == 0)
		return PRECIPITATION;
	else {
		std::cout << "convertDisType #" << dis_type_string << "# not found"
				<< std::endl;
		exit(1);
	}
	return INVALID_DIS_TYPE;
}

std::string convertDisTypeToString(DistributionType dis_type)
{
	if (dis_type == ANALYTICAL) return "ANALYTICAL";
	if (dis_type == AVERAGE) return "AVERAGE";
	if (dis_type == CONSTANT) return "CONSTANT";
	if (dis_type == CONSTANT_GEO) return "CONSTANT_GEO";
	if (dis_type == GRADIENT) return "GRADIENT";
	if (dis_type == RESTART) return "RESTART";
	if (dis_type == LINEAR) return "LINEAR";
	if (dis_type == POINT) return "POINT";
	if (dis_type == CONSTANT_NEUMANN) return "CONSTANT_NEUMANN";
	if (dis_type == LINEAR_NEUMANN) return "LINEAR_NEUMANN";
	if (dis_type == NORMALDEPTH) return "NORMALDEPTH";
	if (dis_type == CRITICALDEPTH) return "CRITICALDEPTH";
	if (dis_type == GREEN_AMPT) return "GREEN_AMPT";
	if (dis_type == SYSTEM_DEPENDENT) return "SYSTEM_DEPENDENT";
	if (dis_type == PRECIPITATION) return "PRECIPITATION";

	return "INVALID_DIS_TYPE";
}
} // end namespace FiniteElement
