/*!
    \file rf_pcs1.cpp
     \brief definitions of member functions of CRFProcess

     for JFNK method
*/

#if defined (NEW_EQS) && defined(JFNK_H2M)

#include <cfloat>
#include <iostream>
#include <iomanip>                                //WW
#include <algorithm>
#include "rf_pcs.h"
#include "pcs_dm.h"
#include "equation_class.h"

using namespace std;
////////////////////////////////////////////////////////////////////////////
///
/// Assemble the residual for JFNK method. 08/2010 WW
///
////////////////////////////////////////////////////////////////////////////

/// Get du from the buffer
void CRFProcess::Recovery_du_JFNK()
{
   int k, p_idx;
   long i, j;

   j = 0;
   for(k=0; k<pcs_number_of_primary_nvals; k++)
   {
      p_idx = p_var_index[k];
      if(k<problem_dimension_dm)
         p_idx--;

      for(i=0; i<num_nodes_p_var[k]; i++)
      {
         SetNodeValue(i, p_idx, array_u_JFNK[j]);
         j++;
      }
   }

   /// Recovery Dirchlet BC
   /// 03.11.2010
   for(i=0; i<(long)BC_JFNK.size(); i++)
   {
      bc_JFNK bc_entry = BC_JFNK[i];
      if(bc_entry.incremental)
         continue;
      SetNodeValue(bc_entry.bc_node, bc_entry.var_idx,
         GetNodeValue(bc_entry.bc_node, bc_entry.var_idx)-bc_entry.bc_value);
   }

}


///////////////////////////////////////////////////////////////////////////////////
///
///
/// Force term control for inexact Newton method.
///
///   01.2011. WW
///////////////////////////////////////////////////////////////////////////////////
bool CRFProcess::ForceTermCriterion(double *Jdx, const int num_iteration)
{
   long i;
   double norm_f_pls_Jdx, f_factor = 0.5;
   bool done = false;

   if(num_iteration)
      f_factor = 0.1;

   /// Calculate ||J*dx+F(x)||
   norm_f_pls_Jdx = 0.;
   double diff_v;
   for(i=0; i<eqs_new->size_global; i++)
   {
      diff_v = Jdx[i] + array_Fu_JFNK[i];
      norm_f_pls_Jdx  += diff_v*diff_v;
   }
   norm_f_pls_Jdx = sqrt(norm_f_pls_Jdx);

   if(norm_f_pls_Jdx < f_factor*norm_F0)
      done = true;

   return done;
}


///////////////////////////////////////////////////////////////////////////////////
///
///
///  Calculate J*v for the JFNK method
///
///   10.2010. WW
///////////////////////////////////////////////////////////////////////////////////
void CRFProcess:: Jacobian_Multi_Vector_JFNK(double *v, double *Jv)
{

   //----------------------------------------------------------------------
   int k, p_idx;
   long i, j;
   //float sign_rhs;
   bool HM = false;
   if(problem_dimension_dm!=pcs_number_of_primary_nvals)
      HM = true;

   /// Calculate F(u): u -->  array_u_JFNK
   if(v == NULL)
   {
      BC_JFNK.clear();

      JFNK_precond = false;

      /// Assemble Jacabin preconditioner
      if(m_num->ls_precond==1)
      {
         JFNK_precond = true;
         eqs_new->Init_Precond_Jacobi_JFNK();
      }

      /// Assemble F(u)--> eqs->f
      IncorporateBoundaryConditions();
      GlobalAssembly();
      for(i=0; i<(long)BC_JFNK.size(); i++)
      {
         bc_JFNK bc_entry = BC_JFNK[i];
         eqs_new->b[bc_entry.bc_eqs_idx] = bc_entry.bc_value;
      }

      j = 0;
      /// Calculate norm_F0 for line search
      norm_F0 = 0.;
      for(k=0; k<pcs_number_of_primary_nvals; k++)
      {
         //sign_rhs = 1.0;
         //if(k>=problem_dimension_dm)
         //   sign_rhs = -1.0;

         for(i=0; i<num_nodes_p_var[k]; i++)
         {
            array_Fu_JFNK[j] = -eqs_new->b[j];
            norm_F0 += array_Fu_JFNK[j]*array_Fu_JFNK[j];
            j++;
         }
      }
      norm_F0 = sqrt(norm_F0);

      norm_u_JFNK[0] = 0.;
      if(HM)
         norm_u_JFNK[1] = 0.;
      j = 0;
      for(k=0; k<pcs_number_of_primary_nvals; k++)
      {
         p_idx = p_var_index[k];
         if(k<problem_dimension_dm&&HM)
            p_idx--;

         for(i=0; i<num_nodes_p_var[k]; i++)
         {
            array_u_JFNK[j] = GetNodeValue(i, p_idx);
            norm_u_JFNK[0] += array_u_JFNK[j]*array_u_JFNK[j];
            j++;
         }
         if((k == problem_dimension_dm-1)&&HM)
            norm_u_JFNK[1] = norm_u_JFNK[0];
      }

      if(HM)
      {
         norm_u_JFNK[0] -= norm_u_JFNK[1];
         norm_u_JFNK[1] = sqrt(norm_u_JFNK[1]);
      }
      norm_u_JFNK[0] = sqrt(norm_u_JFNK[0]);
      return;
   }

   //
   long size;
   double perturbation[2], pert = 0.;
   double norm_v[2];
#define  pert_2
#ifdef pert_2
   double b_fac =  sqrt(DBL_EPSILON);
#endif
   double pert_defalt = sqrt(DBL_EPSILON);
   int num_pcs = 1;

   /// Do not assemble Jacabin preconditioner during linear solver
   JFNK_precond = false;

   norm_v[0] = norm_v[1] = 0.;
   j = 0;
   for(k=0; k<pcs_number_of_primary_nvals; k++)
   {

      for(i=0; i<num_nodes_p_var[k]; i++)
      {
         norm_v[0] += v[j]*v[j];
         j++;
      }
      if((k == problem_dimension_dm-1)&&HM)
         norm_v[1] = norm_v[0];
   }
   if(HM)
   {
      norm_v[0] -= norm_v[1];
      norm_v[1] = sqrt(norm_v[1]);
      num_pcs = 2;
   }
   norm_v[0] = sqrt(norm_v[0]);
   size = j;

   if((norm_v[0]+norm_v[1])>DBL_MIN)
   {

#define apert_a
#ifdef pert_a
      for(k=0; k<num_pcs; k++)
      {
         if(norm_v[k]>DBL_MIN)
         {
#define  apert_2
#ifdef pert_2
            perturbation[k] = 0.;

            for(i=0; i<size; i++)
               perturbation[k] += fabs(array_u_JFNK[i]);
            perturbation[k] = b_fac*perturbation[k]/(norm_v[k]*size)+b_fac;
            if(perturbation[k]<DBL_MIN)
               perturbation[k] = pert_defalt;
#else
            perturbation[k] = sqrt((1.+norm_u_JFNK[k])*sqrt(DBL_EPSILON))/norm_v[k];
#endif
            if(perturbation[k]<DBL_MIN)
               perturbation[k] = pert_defalt;
         }
         else
            perturbation[k] = pert_defalt;
      }
#else                                       //pert_b
      double udv[2], uv[2], suv[2], typu[2], sign = 1.;
      for(k=0; k<num_pcs; k++)
      {
         udv[k] = uv[k] = suv[k] =0.;
         typu[k] =  0.;
      }

      j = 0;
      for(k=0; k<pcs_number_of_primary_nvals; k++)
      {

         for(i=0; i<num_nodes_p_var[k]; i++)
         {
            udv[0] += fabs(array_u_JFNK[j]*v[j]);
            uv[0] += fabs(v[j]);
            suv[0] += fabs(array_u_JFNK[j]*v[j]);

            typu[0] = max(typu[0], fabs(array_u_JFNK[j]));
            j++;
         }
         if((k == problem_dimension_dm-1)&&HM)
         {
            udv[1] = udv[0];
            uv[1] = uv[0];
            suv[1] = suv[0];
            typu[1] = typu[0];
         }
      }
      if(HM)
      {
         udv[0] -= udv[1];
         uv[0] -= uv[1];
         suv[0] -= suv[1];
         swap(typu[0], typu[1]);
      }

      for(k=0; k<num_pcs; k++)
      {
         if(norm_v[k]>DBL_MIN)
         {
            if(fabs(suv[k])>DBL_EPSILON)
               sign = fabs(suv[k])/suv[k];
            else
               sign = 1.0;

            perturbation[k] = sign*b_fac*max(udv[k], typu[k]*uv[k])/norm_v[k];
         }
         else
            perturbation[k] = pert_defalt;
         if(perturbation[k]<DBL_MIN)
            perturbation[k] = pert_defalt;

      }
#endif

      /// Initialize rhs
      for(i=0; i<eqs_new->size_global; i++)
         eqs_new->b[i] = 0.0;                     ///-F(u)

      //TEST
      perturbation[0] = sqrt(DBL_EPSILON);

      ///1. For PDEs excluding that of deformation. 24.22.2010
      j= 0;
      for(k=0; k<pcs_number_of_primary_nvals; k++)
      {
         p_idx = p_var_index[k];
         if(k<problem_dimension_dm&&HM)
            p_idx--;

         for(i=0; i<num_nodes_p_var[k]; i++)
         {
            SetNodeValue(i, p_idx, array_u_JFNK[j]+perturbation[0]*v[j]);
            j++;
         }
      }
      /// Assemble F(u+epsilon*v) and calculate F(u+epsilon*v)-F(u)
      /// Apply Dirchlet BC.  u_b --> node_value
      /// 01.11.2010
      for(i=0; i<(long)BC_JFNK.size(); i++)
      {
         bc_JFNK bc_entry = BC_JFNK[i];
         if(!bc_entry.incremental)
            continue;
         SetNodeValue(bc_entry.bc_node, bc_entry.var_idx, bc_entry.bc_value0);
      }
      GlobalAssembly_std();

      ///2. For the PDE of deformation
      if(HM)
      {
         j= 0;
         for(k=0; k<pcs_number_of_primary_nvals; k++)
         {
            p_idx = p_var_index[k];
            if(k<problem_dimension_dm&&HM)
               p_idx--;

            for(i=0; i<num_nodes_p_var[k]; i++)
            {
               SetNodeValue(i, p_idx, array_u_JFNK[j]+perturbation[1]*v[j]);
               j++;
            }
         }
         /// Assemble F(u+epsilon*v) and calculate F(u+epsilon*v)-F(u)
         /// Apply Dirchlet BC.  u_b --> node_value
         /// 01.11.2010
         for(i=0; i<(long)BC_JFNK.size(); i++)
         {
            bc_JFNK bc_entry = BC_JFNK[i];
            if(bc_entry.incremental)
               continue;
            SetNodeValue(bc_entry.bc_node, bc_entry.var_idx, bc_entry.bc_value);
         }
         GlobalAssembly_DM();
      }

      IncorporateSourceTerms();

      /* 
      /// 24.11.2010. WW
      if(ite_steps>1)
         IncorporateBoundaryConditions();
      */

      j = 0;
      for(k=0; k<pcs_number_of_primary_nvals; k++)
      {
         //sign_rhs = 1.0;
         //if(k>=problem_dimension_dm)
         //   sign_rhs = -1.0;

         if(k<problem_dimension_dm&&HM)
            pert = perturbation[1];
         else
            pert = perturbation[0];

         for(i=0; i<num_nodes_p_var[k]; i++)
         {
                                                  ///Jv
            Jv[j] = (-eqs_new->b[j]-array_Fu_JFNK[j])/pert;
            j++;
         }
      }
      /// Apply Dirchlet BC.  x_i = x_i^0 as
      ///   F(u_u) = F(u_i) x_i = F(u_i)*x_i^0 = b_i
      /// 20.09.2010
      for(i=0; i<(long)BC_JFNK.size(); i++)
      {
         bc_JFNK bc_entry = BC_JFNK[i];
         Jv[bc_entry.bc_eqs_idx] = v[bc_entry.bc_eqs_idx];
      }

   }
   else                                           /// If x is zero.
   {
      j = 0;
      for(k=0; k<pcs_number_of_primary_nvals; k++)
      {

         for(i=0; i<num_nodes_p_var[k]; i++)
         {
            Jv[j] = 0.;
            j++;
         }
      }

   }

}


///////////////////////////////////////////////////////////////////////////////////
///
///
///  Line search for the JFNK method
///
///   10.2010. WW
///////////////////////////////////////////////////////////////////////////////////
double CRFProcess::LineSearch()
{
   int k, p_idx;
   long i, j;
   double damping = 1.0;
   //float sign_rhs;
   bool HM = false;
   if(problem_dimension_dm!=pcs_number_of_primary_nvals)
      HM = true;

   for(;;)
   {
      /// Initialize rhs
      for(i=0; i<eqs_new->size_global; i++)
         eqs_new->b[i] = 0.0;                     ///-F(u)

      j= 0;
      for(k=0; k<pcs_number_of_primary_nvals; k++)
      {
         p_idx = p_var_index[k];
         if(k<problem_dimension_dm&&HM)
            p_idx--;

         for(i=0; i<num_nodes_p_var[k]; i++)
         {
            SetNodeValue(i, p_idx, array_u_JFNK[j]+damping*eqs_new->x[j]);
            j++;
         }
      }
      GlobalAssembly();
      for(i=0; i<(long)BC_JFNK.size(); i++)
      {
         bc_JFNK bc_entry = BC_JFNK[i];
         eqs_new->b[bc_entry.bc_eqs_idx] = bc_entry.bc_value;
      }

      double normFplsX = 0.;
      /// j is the size of the vector
      for(i=0; i<j; i++)
         normFplsX += eqs_new->b[i]*eqs_new->b[i];

      if(sqrt(normFplsX)<norm_F0)
         break;
      else
         damping *= 0.5;
      /// If the damping is too small, jump out
      if(damping<1.e-3)
         break;

   }

   return damping;
}
#endif
