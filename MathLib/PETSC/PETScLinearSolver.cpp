/*!
   \brief Definition of member functions of class PETScLinearSolver

   10~11.2011. WW

*/
#include "PETScLinearSolver.h"

#include <iostream>
#include <list>

#include "../../FEM/display.h"
#include "StringTools.h"

namespace petsc_group
{

 PETScLinearSolver :: PETScLinearSolver (const int size)
   :lsolver(NULL), prec(NULL), global_x0(NULL),  global_x1(NULL), global_buff(NULL) 
 {
   ltolerance = 1.e-10;
   m_size = size;
   time_elapsed = 0.0;   
   d_nz = 10; 
   o_nz = 10; 	
   nz = 10;   
   m_size_loc = PETSC_DECIDE;
   mpi_size = 0;
   rank = 0;
 }

PETScLinearSolver:: ~PETScLinearSolver()
{
  VecDestroy(&b);
  VecDestroy(&x);
  MatDestroy(&A);
  if(lsolver) KSPDestroy(&lsolver);
  // if(prec) PCDestroy(&prec);

  if(global_x0)
    delete []  global_x0;
  if(global_x1)
    delete []  global_x1;
  if(global_buff)
    delete []  global_buff;

  PetscPrintf(PETSC_COMM_WORLD,"\n>>Number of Unknows: %d", m_size);
  PetscPrintf(PETSC_COMM_WORLD,"\n>>Elapsed time in linear solver: %fs", time_elapsed);
}

void PETScLinearSolver::Init(const int *sparse_index)
{
   if(sparse_index)
   {
      d_nz = sparse_index[0]; 
      o_nz = sparse_index[1]; 	
      nz = sparse_index[2]; 
      m_size_loc = sparse_index[3];          
   }   
    
   VectorCreate(m_size);   
   MatrixCreate(m_size, m_size);

   global_x0 = new PetscScalar[m_size];
   global_x1 = new PetscScalar[m_size];
   global_buff = new PetscScalar[m_size];

}


/*!
  \brief KSP and PC type

 KSPRICHARDSON "richardson"
 KSPCHEBYCHEV  "chebychev"
 KSPCG         "cg"
 KSPCGNE       "cgne"
 KSPNASH       "nash"
 KSPSTCG       "stcg"
 KSPGLTR       "gltr"
 KSPGMRES      "gmres"
 KSPFGMRES     "fgmres" 
 KSPLGMRES     "lgmres"
 KSPDGMRES     "dgmres"
 KSPTCQMR      "tcqmr"
 KSPBCGS       "bcgs"
 KSPIBCGS        "ibcgs"
 KSPBCGSL        "bcgsl"
 KSPCGS        "cgs"
 KSPTFQMR      "tfqmr"
 KSPCR         "cr"
 KSPLSQR       "lsqr"
 KSPPREONLY    "preonly"
 KSPQCG        "qcg"
 KSPBICG       "bicg"
 KSPMINRES     "minres"
 KSPSYMMLQ     "symmlq"
 KSPLCD        "lcd"
 KSPPYTHON     "python"
 KSPBROYDEN    "broyden"
 KSPGCR        "gcr"
 KSPNGMRES     "ngmres"
 KSPSPECEST    "specest"

 PCNONE            "none"
 PCJACOBI          "jacobi"
 PCSOR             "sor"
 PCLU              "lu"
 PCSHELL           "shell"
 PCBJACOBI         "bjacobi"
 PCMG              "mg"
 PCEISENSTAT       "eisenstat"
 PCILU             "ilu"
 PCICC             "icc"
 PCASM             "asm"
 PCGASM            "gasm"
 PCKSP             "ksp"
 PCCOMPOSITE       "composite"
 PCREDUNDANT       "redundant"
 PCSPAI            "spai"
 PCNN              "nn"
 PCCHOLESKY        "cholesky"
 PCPBJACOBI        "pbjacobi"
 PCMAT             "mat"
 PCHYPRE           "hypre"
 PCPARMS           "parms"
 PCFIELDSPLIT      "fieldsplit"
 PCTFS             "tfs"
 PCML              "ml"
 PCPROMETHEUS      "prometheus"
 PCGALERKIN        "galerkin"
 PCEXOTIC          "exotic"
 PCHMPI            "hmpi"
 PCSUPPORTGRAPH    "supportgraph"
 PCASA             "asa"
 PCCP              "cp"
 PCBFBT            "bfbt"
 PCLSC             "lsc"
 PCPYTHON          "python"
 PCPFMG            "pfmg"
 PCSYSPFMG         "syspfmg"
 PCREDISTRIBUTE    "redistribute"
 PCSACUSP          "sacusp"
 PCSACUSPPOLY      "sacusppoly"
 PCBICGSTABCUSP    "bicgstabcusp"
 PCSVD             "svd"
 PCAINVCUSP        "ainvcusp"
 PCGAMG            "gamg"

*/
void PETScLinearSolver::Config(const PetscReal tol, const PetscInt maxits, const KSPType lsol, const PCType prec_type, const std::string &misc_setting)
{
   ltolerance = tol;
   sol_type = lsol;
   pc_type = prec_type; 


   KSPCreate(PETSC_COMM_WORLD,&lsolver);
   KSPSetOperators(lsolver, A, A,DIFFERENT_NONZERO_PATTERN);
   KSPSetType(lsolver,lsol);

   KSPGetPC(lsolver, &prec);
   PCSetType(prec, prec_type); //  PCJACOBI); //PCNONE);
   KSPSetTolerances(lsolver,ltolerance, PETSC_DEFAULT, PETSC_DEFAULT, maxits);

   if (!misc_setting.empty()) {
	   PetscPrintf(PETSC_COMM_WORLD, "-> additional PETSc arguments:\n");
	   std::list<std::string> lst = splitString(misc_setting,' ');
	   for (std::list<std::string>::iterator itr=lst.begin(); itr!=lst.end(); ++itr) {
		   // key-value or only key
		   std::string &str1 = *itr;
		   if (str1.find('-')==std::string::npos)
			   continue;
		   ++itr;
		   if (itr==lst.end())
			   break;
		   std::string val = "";
		   std::string &str2 = *itr;
		   if (str2.find('-')==std::string::npos)
				val = str2;
		   else
			   --itr;
		   vec_para.push_back(std::make_pair(str1, val));
		   PetscPrintf(PETSC_COMM_WORLD, "\t %s = %s\n", str1.c_str(), val.c_str());
	   }
   }
   for (std::vector<Para>::iterator itr=vec_para.begin(); itr!=vec_para.end(); ++itr) {
	   PetscOptionsSetValue(itr->first.c_str(),itr->second.c_str());
   }
   char* copts;
   PetscOptionsGetAll(&copts);
   PetscPrintf(PETSC_COMM_WORLD, "-> PETSc options = %s\n", copts);
   PetscFree(copts);
   KSPSetFromOptions(lsolver);
   // reset options
   for (std::vector<Para>::iterator itr=vec_para.begin(); itr!=vec_para.end(); ++itr) {
	   PetscOptionsSetValue(itr->first.c_str(),"");
   }

}
//-----------------------------------------------------------------
void PETScLinearSolver::VectorCreate(PetscInt m)
{
  //PetscErrorCode ierr;  // returned value from PETSc functions 
  VecCreate(PETSC_COMM_WORLD, &b);
  ////VecCreateMPI(PETSC_COMM_WORLD,m_size_loc, m, &b);
  //VecSetSizes(b, m_size_loc, m);
  VecSetSizes(b, PETSC_DECIDE, m);
  VecSetFromOptions(b);
  VecSetUp(b); //kg44 for PETSC 3.3 
  VecDuplicate(b, &x);

  //VecGetOwnershipRange(b, &i_start,&i_end);
}


void PETScLinearSolver::MatrixCreate( PetscInt m, PetscInt n)
{
  PetscErrorCode ierr;
  MatCreate(PETSC_COMM_WORLD, &A);
  // TEST  MatSetSizes(A, m_size_loc, PETSC_DECIDE, m, n);

  ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,m,n);
  //MatSetSizes(A, m_size_loc, PETSC_DECIDE, m,  n);
  CHKERRCONTINUE(ierr);

  MatSetFromOptions(A);
  MatSetType(A,MATMPIAIJ);
  ScreenMessage2("-> set PETSc matrix preallocation wiht d_nz=%d and o_nz=%d\n", d_nz, o_nz);
  MatMPIAIJSetPreallocation(A,d_nz,PETSC_NULL, o_nz,PETSC_NULL);
  //MatSeqAIJSetPreallocation(A,d_nz,PETSC_NULL);
  MatSetOption(A, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE); // for MatZeroRows()
  MatGetOwnershipRange(A,&i_start,&i_end);
  ScreenMessage2("-> PETSc linear solver range: start=%d, end=%d\n", i_start, i_end);

  //  std::cout<<"sub_a  "<<i_start<<";   sub_d "<<i_end<<"\n";
}

void  PETScLinearSolver::getLocalRowColumnSizes(int *m, int *n)
{
  MatGetLocalSize(A, m, n);
}
void  PETScLinearSolver::getOwnerRange(int *start_r, int *end_r)
{
  *start_r = i_start;
  *end_r = i_end;
}

void PETScLinearSolver::Solver()
{
  
   //TEST
#ifdef TEST_MEM_PETSC
   PetscLogDouble mem1, mem2;
   PetscMemoryGetCurrentUsage(&mem1);
#endif
 
  /* 
  //TEST
  PetscViewer viewer;
  PetscViewerASCIIOpen(PETSC_COMM_WORLD, "x.txt", &viewer);
  PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);
  PetscObjectSetName((PetscObject)x,"Solution");
  VecView(x, viewer);   
  */


   int its; 
   PetscLogDouble v1,v2;
   KSPConvergedReason reason;

   // #define PETSC34
   //kg44 quick fix to compile PETSC with version PETSCV3.4
#ifdef USEPETSC34
   PetscTime(&v1);
#else
   PetscGetTime(&v1);
#endif

   PetscPrintf(PETSC_COMM_WORLD,"\n================================================\n");
   PetscPrintf(PETSC_COMM_WORLD, "*** PETSc linear solver\n");
   KSPSolve(lsolver, b, x);
  
   const char *slv_type;
   const char *prc_type;
   KSPGetType(lsolver, &slv_type);
   PCGetType(prec, &prc_type);
   KSPGetConvergedReason(lsolver,&reason); //CHKERRQ(ierr);
   KSPGetIterationNumber(lsolver,&its); //CHKERRQ(ierr);
   PetscReal rtol=.0, abstol=.0, dtol=.0;
   PetscInt maxits = 0;
   KSPGetTolerances(lsolver, &rtol, &abstol, &dtol, &maxits);
   PetscReal r_norm = .0;
   KSPGetResidualNorm(lsolver, &r_norm);
   PetscReal b_norm = .0;
   VecNorm(b, NORM_2, &b_norm);
   PetscReal error_r = r_norm/b_norm;

   PetscPrintf(PETSC_COMM_WORLD, "solver    : %s\n", slv_type);
   PetscPrintf(PETSC_COMM_WORLD, "precon    : %s\n", prc_type);
   PetscPrintf(PETSC_COMM_WORLD, "iteration : %d/%d\n", its, maxits);
   PetscPrintf(PETSC_COMM_WORLD, "residual  : |r|=%e, |b|=%e, error=%e (tol=%e)\n", r_norm, b_norm, error_r, rtol);
   if (reason>=0) {
      PetscPrintf(PETSC_COMM_WORLD, "status    : Converged (reason=%d)\n", reason);
   } else {
      if (reason==KSP_DIVERGED_INDEFINITE_PC)
      {
         PetscPrintf(PETSC_COMM_WORLD, "status    : Diverged (indefinite precon)\n", reason);
         PetscPrintf(PETSC_COMM_WORLD, "            Run the executable again but with -pc_factor_shift_positive_definite option.\n");
      }
      else if (reason==KSP_DIVERGED_ITS)
      {
          PetscPrintf(PETSC_COMM_WORLD, "status    : Diverged (max iteration)\n", reason);
      }
      else
      {
          PetscPrintf(PETSC_COMM_WORLD, "status    : Diverged (reason=%d)\n", reason);
      }
      exit(1);
   }

   PetscPrintf(PETSC_COMM_WORLD,"================================================\n");

   //VecAssemblyBegin(x);
   //VecAssemblyEnd(x);

   //kg44 quick fix to compile PETSC with version PETSCV3.4
#ifdef USEPETSC34
   PetscTime(&v2);
#else
   PetscGetTime(&v2);
#endif

   time_elapsed += v2-v1;

   
#define aTEST_OUT
#ifdef TEST_OUT
  //TEST
   PetscViewer viewer;
   PetscViewerASCIIOpen(PETSC_COMM_WORLD, "x2.txt", &viewer);
   PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);
   PetscObjectSetName((PetscObject)A,"Matrix");
   MatView(A, viewer);
   PetscObjectSetName((PetscObject)x,"Solution");
   VecView(x, viewer);
   PetscObjectSetName((PetscObject)b,"RHS");
   VecView(b, viewer);   
    VecDestroy(&b);
  VecDestroy(&x);
  MatDestroy(&A);
  if(lsolver) KSPDestroy(&lsolver);
  // if(prec) PCDestroy(&prec);
  if(global_x0)
    delete []  global_x0;
  if(global_x1)
    delete []  global_x1;
   PetscFinalize();
   exit(0);
#endif


#ifdef TEST_MEM_PETSC
  //TEST
   PetscMemoryGetCurrentUsage(&mem2);
   PetscPrintf(PETSC_COMM_WORLD, "###Memory usage by solver. Before :%f After:%f Increase:%d\n", mem1, mem2, (int)(mem2 - mem1));
#endif
}

  void PETScLinearSolver::AssembleRHS_PETSc()
{
  VecAssemblyBegin(b);
  VecAssemblyEnd(b);
}
void PETScLinearSolver::AssembleUnkowns_PETSc()
{
  VecAssemblyBegin(x);
  VecAssemblyEnd(x);
}
void PETScLinearSolver::AssembleMatrixPETSc(const MatAssemblyType type)
{
  MatAssemblyBegin(A, type);
  MatAssemblyEnd(A, type);
}


void PETScLinearSolver::UpdateSolutions(PetscScalar *u0, PetscScalar *u1)
{

#ifdef TEST_MEM_PETSC
   //TEST
   PetscLogDouble mem1, mem2;
   PetscMemoryGetCurrentUsage(&mem1);
#endif


  int i, j;
  PetscScalar *xp;
 
  int receivecount;
  PetscInt low,high,otherlow;
  MPI_Status status; 
  PetscInt count;
  int tag = 9999;
  VecGetOwnershipRange(x, &low, &high);
  VecGetLocalSize(x, &count);


  VecGetArray(x, &xp);
  for(i=0; i<count; i++)
    u1[i] = xp[i];


  //double *global_buff = new double[m_size];


  // Collect solution from processes.
  for(j=0; j<count; j++)
    global_buff[low+j] = u1[j];
  for(i=0;i<mpi_size;i++) 
  {
     if(i != rank)
     {

       MPI_Sendrecv( &count, 1, MPI_INT, i,tag, 
                     &receivecount,1,MPI_INT,i,tag, PETSC_COMM_WORLD ,&status);
       MPI_Sendrecv( &low, 1, MPI_INT, i,tag,
                 &otherlow,1,MPI_INT,i,tag,PETSC_COMM_WORLD,&status );
       MPI_Sendrecv( u1, count, MPI_DOUBLE, i,tag,
                     u0,receivecount,MPI_DOUBLE,i,tag, PETSC_COMM_WORLD,&status  );
       for(j=0;j<receivecount;j++)
         global_buff[otherlow+j] = u0[j];
     }
  }


  //MPI_Barrier(PETSC_COMM_WORLD);
  // Copy the collected solution to the array for the previous solution
  for(i=0;i<m_size;i++)
  {
    u1[i] = global_buff[i];
    u0[i] = global_buff[i];
  }
 

  //delete [] global_buff;

  VecRestoreArray(x, &xp);


  //TEST
#ifdef TEST_MEM_PETSC
   PetscMemoryGetCurrentUsage(&mem2);
   PetscPrintf(PETSC_COMM_WORLD, "### Memory usage by Updating. Before :%f After:%f Increase:%d\n", mem1, mem2, (int)(mem2 - mem1));
#endif

}

void PETScLinearSolver::MappingSolution()
{
  UpdateSolutions(global_x0, global_x1);
}


int PETScLinearSolver::GetLocalSolution(PetscScalar *x_l)
{
  PetscInt count;
  VecGetLocalSize(x, &count);

  VecGetArray(x, &x_l);

  return count;
}


int PETScLinearSolver::GetLocalRHS(PetscScalar *rhs_l)
{
  PetscInt count;
  VecGetLocalSize(b, &count);

  VecGetArray(b, &rhs_l);

  return count;
}

double *PETScLinearSolver::GetGlobalSolution() const
{
  return global_x1;
}

/*!
  Get values of the specified elements from a global vector

  @param v_type - Indicator for vector: 0: x; 1: rhs
  @param ni 	- number of elements to get
  @param ix 	- indices where to get them from (in global 1d numbering) 
*/
void  PETScLinearSolver::GetVecValues(const int v_type, PetscInt ni,const PetscInt ix[], 
				      PetscScalar y[]) const
{
  if(v_type == 0)
    VecGetValues(x, ni, ix, y);
  else 
    VecGetValues(b, ni, ix, y);
}

/*!
    Get norm of RHS
    @param nmtype  - norm type
                     NORM_1 denotes sum_i |x_i|
                     NORM_2 denotes sqrt(sum_i (x_i)^2)
                     NORM_INFINITY denotes max_i |x_i| 
    06.2012. WW
*/
PetscReal PETScLinearSolver::GetVecNormRHS(NormType  nmtype)
{
  PetscReal norm = 0.;
  VecNorm(b, nmtype, &norm); 
  return norm; 
}
/*!
    Get norm of x
    @param nmtype  - norm type
                     NORM_1 denotes sum_i |x_i|
                     NORM_2 denotes sqrt(sum_i (x_i)^2)
                     NORM_INFINITY denotes max_i |x_i| 
    06.2012. WW
*/
PetscReal PETScLinearSolver::GetVecNormX(NormType  nmtype)
{
  PetscReal norm = 0.;
  VecNorm(x, nmtype, &norm); 
  return norm; 
}


void  PETScLinearSolver::RestoreLocalSolutionArray(PetscScalar *x_l)
{
   VecRestoreArray(x, &x_l);
}
void  PETScLinearSolver::RestoreLocalRHSArray(PetscScalar *rhs_l)
{
   VecRestoreArray(b, &rhs_l);
}

void PETScLinearSolver::set_bVectorEntry(const int i, const double value )
{

  VecSetValues(b,1,&i,&value,INSERT_VALUES);
}
void PETScLinearSolver::set_xVectorEntry(const int i, const double value)
{

  VecSetValues(x,1,&i,&value,INSERT_VALUES);
}

void  PETScLinearSolver::setArrayValues(int arr_idx, PetscInt ni, const PetscInt ix[], 
                                       const PetscScalar y[],InsertMode iora) 
{
   if(arr_idx == 0)
     VecSetValues(x, ni, ix, y, iora); 
   else if(arr_idx == 1)
     VecSetValues(b, ni, ix, y, iora); 
}



void PETScLinearSolver::add_bVectorEntry(const int i, const double value,InsertMode mode )
{

  VecSetValue(b, i, value, mode);
}
void PETScLinearSolver::add_xVectorEntry(const int i, const double value, InsertMode mode)
{

  VecSetValue(x, i, value,mode);
}


void PETScLinearSolver::Initialize( )
{

   VecSet(b, 0.0);
   VecSet(x, 0.0);
   MatZeroEntries(A);
} 



void PETScLinearSolver::addMatrixEntry(const int i, const int j, const double value)
{

  MatSetValue(A, i, j, value, ADD_VALUES);
}

void PETScLinearSolver::addMatrixEntries(const int m,const int idxm[], const int n, 
             const int idxn[],const PetscScalar v[])
{
  MatSetValues(A, m, idxm, n, idxn, v, ADD_VALUES);
}




void PETScLinearSolver::zeroRows_in_Matrix(const int nrows, const  PetscInt *rows)
{
  PetscScalar one = 1.0;
  // Each process indicates only rows it owns that are to be zeroed
  // MatSetOption(A, MAT_NO_OFF_PROC_ZERO_ROWS,PETSC_TRUE);
  if(nrows>0)
    MatZeroRows (A, nrows, rows, one, PETSC_NULL, PETSC_NULL);
  else
    MatZeroRows(A, 0, PETSC_NULL, one, PETSC_NULL, PETSC_NULL);
}



void PETScLinearSolver::EQSV_Viewer(std::string file_name)
{
  PetscViewer viewer;
  std::string fname = file_name + "_eqs_dump.txt";
  PetscViewerASCIIOpen(PETSC_COMM_WORLD, fname.c_str(), &viewer);
  PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);

  AssembleRHS_PETSc();
  AssembleUnkowns_PETSc();
  AssembleMatrixPETSc(MAT_FINAL_ASSEMBLY );


  //PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_VTK);
  PetscObjectSetName((PetscObject)A,"Stiffness_matrix");
  PetscObjectSetName((PetscObject)b,"RHS");
  PetscObjectSetName((PetscObject)x,"Solution");
  MatView(A,viewer);
  VecView(b, viewer);
  VecView(x, viewer);  

//#define  EXIT_TEST
#ifdef EXIT_TEST 
  VecDestroy(&b);
  VecDestroy(&x);
  MatDestroy(&A);
  if(lsolver) KSPDestroy(&lsolver);
  // if(prec) PCDestroy(&prec);
  if(global_x0)
    delete []  global_x0;
  if(global_x1)
    delete []  global_x1;
   PetscFinalize();
   exit(0);
#endif 
 
}

} //end of namespace
