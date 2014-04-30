
#include "rf_pcs_TH.h"

#ifdef USE_PETSC
#include <petscksp.h>
#endif

#include "StringTools.h"
#if defined(NEW_EQS)
#include "equation_class.h"
#elif defined(USE_PETSC)
#include "PETSC/PETScLinearSolver.h"
#endif

#include "fem_ele_std.h"

CRFProcessTH::CRFProcessTH()
	: CRFProcess(), error_k0(0.0)
{}

CRFProcessTH::~CRFProcessTH()
{
}

void CRFProcessTH::Initialization()
{
	const size_t bufferSize = GetPrimaryVNumber() * m_msh->GetNodesNumber(false);
	if(m_num->nls_method != FiniteElement::NL_JFNK)
		ARRAY.resize(bufferSize);
	int Axisymm = 1; // ani-axisymmetry
	if (m_msh->isAxisymmetry())
		Axisymm = -1;  // Axisymmetry is true
	fem = new CFiniteElementStd(this, Axisymm * m_msh->GetCoordinateFlag());
	fem->SetGaussPointNumber(m_num->ele_gauss_points);

#if 1
	// calculate initial velocity
	ScreenMessage("-> compute velocity from initial pressure\n");
	CalIntegrationPointValue();
#endif
}

#ifdef USE_PETSC
void CRFProcessTH::setSolver( petsc_group::PETScLinearSolver *petsc_solver )
{
   eqs_new = petsc_solver;

	if (this->m_num->petsc_split_fields) {
		ScreenMessage("-> prepare field splits in PETSc\n");
		PetscErrorCode ierr;
		int dof = this->GetPrimaryVNumber();
		eqs_new->vec_subA.resize(dof*dof);
		eqs_new->vec_isg.resize(dof);
		eqs_new->vec_subRHS.resize(dof);
		// setup matrix
		int n_nodes = this->m_msh->getNumNodesGlobal();
		for (int i=0; i<dof*dof; i++) {
			ierr = MatCreate(PETSC_COMM_WORLD,&eqs_new->vec_subA[i]);CHKERRCONTINUE(ierr);
			std::string str = "a" + number2str(i/dof) + number2str(i%dof) + "_";
			ierr = MatSetOptionsPrefix(eqs_new->vec_subA[i],str.c_str());CHKERRCONTINUE(ierr);
			ierr = MatSetSizes(eqs_new->vec_subA[i], PETSC_DECIDE,PETSC_DECIDE,n_nodes,n_nodes);CHKERRCONTINUE(ierr);
			ierr = MatSetType(eqs_new->vec_subA[i], MATMPIAIJ);CHKERRCONTINUE(ierr);
			ierr = MatSetFromOptions(eqs_new->vec_subA[i]);CHKERRCONTINUE(ierr);
			ierr = MatMPIAIJSetPreallocation(eqs_new->vec_subA[i], eqs_new->d_nz, PETSC_NULL, eqs_new->o_nz, PETSC_NULL);CHKERRCONTINUE(ierr);
			MatSetOption(eqs_new->vec_subA[i],MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE);
			MatSetOption(eqs_new->vec_subA[i], MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE); // for MatZeroRows()
			int i_start, i_end;
			MatGetOwnershipRange(eqs_new->vec_subA[i],&i_start,&i_end);
//			ScreenMessage2("-> Sub A[%d]: start=%d, end=%d\n", i, i_start, i_end);
		}
		ierr = MatDestroy(&eqs_new->A);CHKERRCONTINUE(ierr);
		ierr = VecDestroy(&eqs_new->b);CHKERRCONTINUE(ierr);
		ierr = VecDestroy(&eqs_new->x);CHKERRCONTINUE(ierr);
		ierr = MatCreateNest(PETSC_COMM_WORLD, dof, NULL, dof, NULL, &eqs_new->vec_subA[0], &eqs_new->A);CHKERRCONTINUE(ierr);
		ierr = MatGetVecs(eqs_new->A,&eqs_new->b,&eqs_new->x);CHKERRCONTINUE(ierr);
		ierr = VecDuplicate(eqs_new->x, &eqs_new->total_x);CHKERRCONTINUE(ierr);
		// setup index sets
		ierr = MatNestGetISs(eqs_new->A, &eqs_new->vec_isg[0], NULL);CHKERRCONTINUE(ierr);
#if 0
		for (int i=0; i<dof; i++) {
			int is_global_size, is_local_size;
			ISGetSize(eqs_new->vec_isg[i], &is_global_size);
			ISGetLocalSize(eqs_new->vec_isg[i], &is_local_size);
			ScreenMessage2("-> IS[%d]: global size=%d, local size=%d\n", i, is_global_size, is_local_size);
		}
#endif
//		ISView(eqs_new->vec_isg[0],PETSC_VIEWER_STDOUT_WORLD);CHKERRCONTINUE(ierr);
//		ISView(eqs_new->vec_isg[1],PETSC_VIEWER_STDOUT_WORLD);CHKERRCONTINUE(ierr);
	}

	eqs_new->Config(m_num->ls_error_tolerance, m_num->ls_max_iterations, m_num->getLinearSolverName(), m_num->getPreconditionerName(), m_num->ls_extra_arg);

//	if (this->m_num->petsc_split_fields) {
//		PetscErrorCode ierr;
//		int dof = this->GetPrimaryVNumber();
//		// setup pc
//		for (int i=0; i<dof; i++) {
//			ierr = PCFieldSplitSetIS(eqs_new->prec, number2str(i).c_str(), eqs_new->vec_isg[i]);CHKERRCONTINUE(ierr);
//		}
//	}
}
#endif

bool CRFProcessTH::checkNRConvergence()
{
	return true;
#if 0
#ifdef USE_PETSC
	const double g_nnodes = m_msh->getNumNodesLocal();
#else
	const double g_nnodes = m_msh->GetNodesNumber(false);
#endif
	std::vector<double> pval_error(pcs_number_of_primary_nvals);
	unsigned num_dof_errors = 0;

	switch(m_num->getNonLinearErrorMethod())
	{
		// --> ERNORM:	|x1-x0|/|x0|
		//     Norm of the solution vector delta divided by the solution vector (relative error).
		//     A single tolerance applied to all primary variables.
		//
		case FiniteElement::ERNORM:
		{
#if defined(NEW_EQS)
			const double NormDx = eqs_new->NormX();
#elif !defined(USE_PETSC)
			const double NormDx = NormOfUnkonwn_orRHS();
#endif

			double NormX0 = 0.0;
			for(int ii=0;ii<pcs_number_of_primary_nvals;ii++)
			{
				int nidx1 = GetNodeValueIndex(pcs_primary_function_name[ii]) + 1;
				//
				for(long i = 0; i < g_nnodes; i++){
					long k = m_msh->Eqs2Global_NodeIndex[i];
					double val2 = GetNodeValue(k, nidx1);
					NormX0 += val2*val2;
				}
			}
#ifdef USE_PETSC
			double local_value = NormX0;
			MPI_Allreduce(&NormX0, &local_value, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif
			num_dof_errors = 1;
			pval_error[0] = NormDx / (sqrt(NormX0)+std::numeric_limits<double>::epsilon());
			break;
		}
		//
		// --> LMAX:	max(x1-x0)
		//     Local max error (across all elements) of solution vector delta (absolute error).
		//     Tolerance required for each primary variable.
		//
		case FiniteElement::LMAX:
		{
			double NormDx = .0;
			for(int ii=0;ii<pcs_number_of_primary_nvals;ii++)
			{
				NormDx = 0.0;
				for (long i = 0; i < g_nnodes; i++) {
					NormDx = std::max(NormDx, fabs(eqs_x[i+ii*g_nnodes]));
				}
				pval_error[ii] = NormDx;
			}
			break;
		}
		default:
			break;
	}


	// Norm of residual (RHS in eqs)
#ifdef USE_MPI
	const double NormR = dom->eqsH->NormRHS();
#else
#if defined(NEW_EQS)
	const double NormR = eqs_new->NormRHS();
#elif !defined(USE_PETSC)
	const double NormR = NormOfUnkonwn_orRHS(false);
#endif
#endif

	// Norm of dx
#if defined(NEW_EQS)
	const double NormDx = eqs_new->NormX();
#elif !defined(USE_PETSC)
	const double NormDx = NormOfUnkonwn_orRHS();
#endif

	// Check the convergence
	Error1 = Error;
	ErrorU1 = ErrorU;
	if(ite_steps == 1 && this->first_coupling_iteration)
	{
		InitialNorm = NormR;
		InitialNormU0 = NormDx;
		static bool firstime = true;
		if (firstime) {
			InitialNormU = NormDx;
			firstime = false;
		}
	}

	Error = NormR / InitialNorm;
	ErrorU = NormDx / InitialNormU0;
	if(NormR < newton_tol && Error > NormR)
		Error = NormR;
	//           if(Norm<TolNorm)  Error = 0.01*Tolerance_global_Newton;
	if((NormDx / InitialNormU) <= newton_tol)
		Error = NormDx / InitialNormU;

	// Compute damping for Newton-Raphson step
	double damping = 1.0;
	//           if(Error/Error1>1.0e-1) damping=0.5;
	if(Error / Error1 > 1.0e-1 || ErrorU / ErrorU1 > 1.0e-1)
		damping = 0.5;
	if(ErrorU < Error)
		Error = ErrorU;
	// JT: Store the process and coupling errors
	pcs_num_dof_errors = 1;
	if (ite_steps == 1){
		pcs_absolute_error[0] = NormDx;
		pcs_relative_error[0] = pcs_absolute_error[0] / newton_tol;
		cpl_max_relative_error = pcs_relative_error[0];
		cpl_num_dof_errors = 1;
	} else {
		pcs_absolute_error[0] = Error;
		pcs_relative_error[0] = Error / newton_tol;
	}
	//
#ifdef USE_MPI
	if(myrank == 0)
	{
#endif
	//Screan printing:
	std::cout<<"      -->End of Newton-Raphson iteration: "<<ite_steps<<"/"<< n_max_iterations <<std::endl;
	std::cout.width(8);
	std::cout.precision(2);
	std::cout.setf(std::ios::scientific);
	std::cout<<"         NR-Error"<<"  "<<"RHS Norm 0"<<"  "<<"RHS Norm  "<<"  "<<"Unknowns Norm"<<"  "<<"Damping"<<std::endl;
	std::cout<<"         "<<Error<<"  "<<InitialNorm<<"  "<<NormR<<"   "<<NormDx<<"   "<<"   "<<damping<<std::endl;
	std::cout <<"      ------------------------------------------------"<<std::endl;
#ifdef USE_MPI
}
#endif
	if(Error > 100.0 && ite_steps > 1)
	{
		printf ("\n  Attention: Newton-Raphson step is diverged. Programme halt!\n");
		exit(1);
	}
	if(InitialNorm < 10 * newton_tol
		|| NormR < 0.001 * InitialNorm
		|| Error <= newton_tol)
		break;
#endif
}


/*************************************************************************
   ROCKFLOW - Function: CRFProcess::
   Task:  Solve plastic deformation by generalized Newton-Raphson method
   Programming:
   02/2003 OK Implementation
   05/2003 WW Polymorphism function by OK
   last modified: 23.05.2003
 **************************************************************************/
double CRFProcessTH::Execute(int loop_process_number)
{
	ScreenMessage("\n================================================\n");
	ScreenMessage("->Process %d : %s\n", loop_process_number, convertProcessTypeToString (getProcessType()).c_str());
	ScreenMessage("================================================\n");

	clock_t dm_time = -clock();

	m_msh->SwitchOnQuadraticNodes(false);
	if(hasAnyProcessDeactivatedSubdomains || NumDeactivated_SubDomains > 0)
		CheckMarkedElement();

	// system matrix
#if defined(USE_PETSC)
	eqs_x = eqs_new->GetGlobalSolution();
	{
		static bool first = true;
		if (first) {
			if (m_num->petsc_split_fields) {
				// set initial value
				Vec sub_x;
				const long number_of_nodes = num_nodes_p_var[0];
				std::vector<int> vec_pos(number_of_nodes);
				for (long j = 0; j < number_of_nodes; j++)
					vec_pos[j] = m_msh->Eqs2Global_NodeIndex[j];
				std::vector<double> vec_values(number_of_nodes);
				const int ColIndex_p = p_var_index[0];
				VecGetSubVector(eqs_new->total_x, eqs_new->vec_isg[0], &sub_x);
				for (long j = 0; j < number_of_nodes; j++)
					vec_values[j] = GetNodeValue(j,ColIndex_p);
				VecSetValues(eqs_new->total_x, number_of_nodes, &vec_pos[0], &vec_values[0], INSERT_VALUES);
				VecAssemblyBegin(sub_x);
				VecAssemblyEnd(sub_x);
				const int ColIndex_T = p_var_index[1];
				VecRestoreSubVector(eqs_new->total_x, eqs_new->vec_isg[0], &sub_x);
				VecGetSubVector(eqs_new->total_x, eqs_new->vec_isg[1], &sub_x);
				for (long j = 0; j < number_of_nodes; j++)
					vec_values[j] = GetNodeValue(j,ColIndex_T);
				VecAssemblyBegin(sub_x);
				VecAssemblyEnd(sub_x);
				VecRestoreSubVector(eqs_new->total_x, eqs_new->vec_isg[1], &sub_x);
				VecAssemblyBegin(eqs_new->total_x);
				VecAssemblyEnd(eqs_new->total_x);
			}
		}
	}
#endif
#if defined(NEW_EQS)                              //WW
	//
#if defined(USE_MPI)
	CPARDomain* dom = dom_vector[myrank];
	long global_eqs_dim = pcs_number_of_primary_nvals * m_msh->GetNodesNumber(true);
	dom->ConfigEQS(m_num, global_eqs_dim, true);
#else
	eqs_new->ConfigNumerics(m_num);       //27.11.2007 WW
#endif
	//
#elif !defined (USE_PETSC)
	SetZeroLinearSolver(eqs);
#endif

#if 0
	if(this->first_coupling_iteration && m_num->nls_method != FiniteElement::NL_JFNK)
		StoreLastSolution();      //u_n-->temp
	ResetCouplingStep();

	// Initialize incremental displacement: w=0
	InitializeNewtonSteps();
#endif


	// Begin Newton-Raphson steps
	double Error = 1.0;
//	double Error1 = 0.0;
	double ErrorU = 1.0;
//	double ErrorU1 = 0.0;
	double InitialNorm = 0.0;
	double InitialNormDx = 0.0;
//	double InitialNormU = 0.0;
	static double rp0 = .0, rT0 = 0;
	static double rp0_L2 = .0, rT0_L2 = 0;
	double dp_max = std::numeric_limits<double>::max(), dT_max = std::numeric_limits<double>::max();
	double dp_L2 = std::numeric_limits<double>::max(), dT_L2 = std::numeric_limits<double>::max();
	double p_max = std::numeric_limits<double>::max(), T_max = std::numeric_limits<double>::max();
	double p_L2 = std::numeric_limits<double>::max(), T_L2 = std::numeric_limits<double>::max();
	double NormDx = std::numeric_limits<double>::max();

	const double newton_tol = m_num->nls_error_tolerance[0];
	const double tol_dp = m_num->nls_error_tolerance[1];
	const double tol_dT = m_num->nls_error_tolerance[2];
	const int n_max_iterations = m_num->nls_max_iterations;

	iter_nlin = 0;
	bool converged = false;
	while(iter_nlin < n_max_iterations)
	{
		iter_nlin++;

		ScreenMessage("------------------------------------------------\n");
		ScreenMessage("-> Nonlinear iteration: %d/%d\n", iter_nlin-1, n_max_iterations);
		ScreenMessage("------------------------------------------------\n");

		//----------------------------------------------------------------------
		// Solve
		//----------------------------------------------------------------------
		// Refresh solver
#if defined(NEW_EQS)
#ifndef USE_MPI
		eqs_new->Initialize(); //27.11.2007 WW
#endif
#elif defined(USE_PETSC) // || defined(other parallel libs)//03.3012. WW
		eqs_new->Initialize();
#else
		SetZeroLinearSolver(eqs);
#endif

		ScreenMessage("-> Assembling equation system...\n");
		GlobalAssembly();

		//
#ifdef USE_MPI
		const double NormR = dom->eqsH->NormRHS();
#elif defined(NEW_EQS)
		const double NormR = eqs_new->NormRHS();
#elif defined(USE_PETSC)
		const double NormR = eqs_new->GetVecNormRHS();
#else
		const double NormR = NormOfUnkonwn_orRHS(false);
#endif
		double rp_max = std::numeric_limits<double>::max(), rT_max = std::numeric_limits<double>::max();
		double rp_L2 = std::numeric_limits<double>::max(), rT_L2 = std::numeric_limits<double>::max();
#if defined(USE_PETSC)
		if (m_num->petsc_split_fields) {
			Vec sub_x;
			VecGetSubVector(eqs_new->b, eqs_new->vec_isg[0], &sub_x);
			VecNorm(sub_x, NORM_2, &rp_L2);
			VecNormBegin(sub_x, NORM_INFINITY, &rp_max);
			VecNormEnd(sub_x, NORM_INFINITY, &rp_max);
			VecRestoreSubVector(eqs_new->b, eqs_new->vec_isg[0], &sub_x);
			VecGetSubVector(eqs_new->b, eqs_new->vec_isg[1], &sub_x);
			VecNorm(sub_x, NORM_2, &rT_L2);
			VecNormBegin(sub_x, NORM_INFINITY, &rT_max);
			VecNormEnd(sub_x, NORM_INFINITY, &rT_max);
			VecRestoreSubVector(eqs_new->b, eqs_new->vec_isg[1], &sub_x);
		}
#endif
		if(iter_nlin == 1 && this->first_coupling_iteration)
		{
			InitialNorm = NormR;
			rp0 = rp_max;
			rT0 = rT_max;
			rp0_L2 = rp_L2;
			rT0_L2 = rT_L2;
			static bool firstime = true;
			if (firstime) {
				firstime = false;
			}
		}
		Error = std::max(rp_max / rp0, rT_max / rT0); //NormR / InitialNorm;
		const double Error_L2 = std::max(rp_L2 / rp0_L2, rT_L2 / rT0_L2);
		const double dx_i = std::max(dp_max/p_max, dT_max/T_max);
		const double dx_L2 = std::max(dp_L2/p_L2, dT_L2/T_L2);
		ScreenMessage("-> Newton-Raphson: r_i=%.3e, r_2=%.3e, dx_i=%.3e, dx_2=%.3e (tol=%g) %d/%d \n", Error, Error_L2, dx_i, dx_L2, newton_tol, iter_nlin-1, n_max_iterations);
		ScreenMessage("|r0|=%.3e, |r|=%.3e, |r|/|r0|=%.3e, |dx|=%.3e\n", InitialNorm, NormR, NormR / InitialNorm, NormDx);
		ScreenMessage("|rp|_i=%.3e, |rT|_i=%.3e, |rp/r0|_i=%.3e, |rT/r0|_i=%.3e\n", rp_max, rT_max, rp_max / rp0, rT_max / rT0);
		ScreenMessage("|rp|_2=%.3e, |rT|_2=%.3e, |rp/r0|_2=%.3e, |rT/r0|_2=%.3e\n", rp_L2, rT_L2, rp_L2 / rp0_L2, rT_L2 / rT0_L2);
		ScreenMessage("|dp|_i=%.3e, |dT|_i=%.3e, |dp/p|_i=%.3e, |dT/T|_i=%.3e (tol.p=%.1e,T=%.1e)\n", dp_max, dT_max, dp_max/p_max, dT_max/T_max, tol_dp, tol_dT);
		ScreenMessage("|dp|_2=%.3e, |dT|_2=%.3e, |dp/p|_2=%.3e, |dT/T|_2=%.3e (tol.p=%.1e,T=%.1e)\n", dp_L2, dT_L2, dp_L2/p_L2, dT_L2/T_L2, tol_dp, tol_dT);
		if (Error < newton_tol) {
			ScreenMessage("-> Newton-Raphson converged\n");
			converged = true;
			break;
		}

		ScreenMessage("-> Calling linear solver...\n");
		// Linear solver
#if defined(USE_MPI)
		dom->eqsH->Solver(eqs_new->x, global_eqs_dim);
#elif defined(NEW_EQS)
#ifdef LIS
		bool compress_eqs = (type/10==4 || this->NumDeactivated_SubDomains>0);
		iter_lin = eqs_new->Solver(this->m_num, compress_eqs); //NW
#else
		iter_lin = eqs_new->Solver(); //27.11.2007
#endif
#elif defined(USE_PETSC)
//		if (write_leqs) {
//			std::string fname = FileName + "_" + convertProcessTypeToString(this->getProcessType()) + "_leqs_assembly.txt";
//			eqs_new->EQSV_Viewer(fname);
//		}
		iter_lin = eqs_new->Solver();
//		if (write_leqs) {
//			std::string fname = FileName + "_" + convertProcessTypeToString(this->getProcessType()) + "_leqs_solution.txt";
//			eqs_new->EQSV_Viewer(fname);
//		}
//		if (iter_nlin==1) {
//			std::string fname = FileName + "_" + convertProcessTypeToString(this->getProcessType()) + "_leqs_residual.txt";
//			eqs_new->Residual_Viewer(fname);
//		}
		if (!m_num->petsc_split_fields) {
			eqs_new->MappingSolution();
			eqs_x = eqs_new->GetGlobalSolution();
		} else {
			VecAXPY(eqs_new->total_x, 1.0, eqs_new->x);
		}
#else
		iter_lin = ExecuteLinearSolver();
#endif
		if (iter_lin<0) {
			accepted = false;
			Tim->last_dt_accepted = false;
			break;
		}
		iter_lin_max = std::max(iter_lin_max, iter_lin);

		//----------------------------------------------------------------------
		// Check convergence
		//----------------------------------------------------------------------
		if (m_num->petsc_split_fields) {
#if defined(USE_PETSC)
			Vec sub_x;
			// dx
			VecGetSubVector(eqs_new->x, eqs_new->vec_isg[0], &sub_x);
			VecNorm(sub_x, NORM_2, &dp_L2);
			VecNormBegin(sub_x, NORM_INFINITY, &dp_max);
			VecNormEnd(sub_x, NORM_INFINITY, &dp_max);
			VecRestoreSubVector(eqs_new->x, eqs_new->vec_isg[0], &sub_x);
			VecGetSubVector(eqs_new->x, eqs_new->vec_isg[1], &sub_x);
			VecNorm(sub_x, NORM_2, &dT_L2);
			VecNormBegin(sub_x, NORM_INFINITY, &dT_max);
			VecNormEnd(sub_x, NORM_INFINITY, &dT_max);
			VecRestoreSubVector(eqs_new->x, eqs_new->vec_isg[1], &sub_x);
			// total x
			VecGetSubVector(eqs_new->total_x, eqs_new->vec_isg[0], &sub_x);
			VecNorm(sub_x, NORM_2, &p_L2);
			VecNormBegin(sub_x, NORM_INFINITY, &p_max);
			VecNormEnd(sub_x, NORM_INFINITY, &p_max);
			VecRestoreSubVector(eqs_new->total_x, eqs_new->vec_isg[0], &sub_x);
			VecGetSubVector(eqs_new->total_x, eqs_new->vec_isg[1], &sub_x);
			VecNorm(sub_x, NORM_2, &T_L2);
			VecNormBegin(sub_x, NORM_INFINITY, &T_max);
			VecNormEnd(sub_x, NORM_INFINITY, &T_max);
			VecRestoreSubVector(eqs_new->total_x, eqs_new->vec_isg[1], &sub_x);
			// relative error
#endif
		}

		// Norm of dx
#if defined(NEW_EQS)
		NormDx = eqs_new->NormX();
#elif defined(USE_PETSC)
		NormDx = eqs_new->GetVecNormX();
#else
		NormDx = NormOfUnkonwn_orRHS();
#endif

		// Check the convergence
//		Error1 = Error;
//		ErrorU1 = ErrorU;
		if(iter_nlin == 1 && this->first_coupling_iteration)
		{
			InitialNormDx = NormDx;
			static bool firstime = true;
			if (firstime) {
//				InitialNormU = NormDx;
				firstime = false;
			}
		}

		ErrorU = NormDx / InitialNormDx;
#if 0
		if(NormR < newton_tol && Error > NormR)
			Error = NormR;
		//           if(Norm<TolNorm)  Error = 0.01*Tolerance_global_Newton;
		if((NormDx / InitialNormU) <= newton_tol)
			Error = NormDx / InitialNormU;
		if(ErrorU < Error)
			Error = ErrorU;
#endif
		// JT: Store the process and coupling errors
		pcs_num_dof_errors = 1;
		if (iter_nlin == 1){
			pcs_absolute_error[0] = NormDx;
			pcs_relative_error[0] = pcs_absolute_error[0] / newton_tol;
			cpl_max_relative_error = pcs_relative_error[0];
			cpl_num_dof_errors = 1;
		} else {
			pcs_absolute_error[0] = Error;
			pcs_relative_error[0] = Error / newton_tol;
		}
		//
		//Screan printing:
//		ScreenMessage("-> update solutions\n");
#if 0
		if(Error > 100.0 && iter_nlin > 1)
		{
			ScreenMessage ("\n  Attention: Newton-Raphson step is diverged. Programme halt!\n");
			accepted = false;
			Tim->last_dt_accepted = false;
			return -1;
		}
#endif
		if (std::max(rp_max / rp0, rT_max / rT0) < newton_tol || (dp_max < tol_dp && dT_max < tol_dT) ) {
			ScreenMessage("-> Newton-Raphson converged\n");
			converged = true;
			break;
		}
//		if(InitialNorm < 10 * newton_tol
//			|| NormR < 0.001 * InitialNorm
//			|| Error <= newton_tol)
//			break;

		// x^k1 = x^k + dx
		UpdateIterativeStep(1.0);

		ScreenMessage("-> update velocity\n");
		CalIntegrationPointValue();
	} // Newton-Raphson iteration

	// x^k1 = x^k + dx
	UpdateIterativeStep(1.0);

	iter_nlin_max = std::max(iter_nlin_max, iter_nlin);

	if (!converged && m_num->nls_max_iterations>1)
		accepted = false;

	//
	dm_time += clock();
#if defined(USE_MPI) || defined(USE_PETSC)
	if(myrank == 0)
	{
#endif
    std::cout <<"CPU time elapsed in this process: " << (double)dm_time / CLOCKS_PER_SEC<<"s"<<std::endl;
    std::cout <<"------------------------------------------------"<<std::endl;
#if defined(USE_MPI) || defined(USE_PETSC)
}
#endif
	// Recovery the old solution.  Temp --> u_n	for flow proccess
//	RecoverSolution();
	//
#ifdef NEW_EQS                              //WW
#if defined(USE_MPI)
	dom->eqsH->Clean();
#else
	// Also allocate temporary memory for linear solver. WW
	eqs_new->Clean();
#endif
#endif

	// For coupling control
	Error = 0.0;
	for(size_t n = 0; n < m_msh->GetNodesNumber(false); n++)
	{
		for(int l = 0; l < pcs_number_of_primary_nvals; l++)
		{
			double NormU = GetNodeValue(n, 2*l) - GetNodeValue(n, 2*l+1);
			Error += NormU * NormU;
		}
	}
	double sqrt_norm = sqrt(Error);

	if(this->first_coupling_iteration)
		Error = fabs(sqrt_norm - error_k0) / error_k0;
	else
		Error = sqrt_norm;
	error_k0 = sqrt_norm;

	if(!accepted || Tim->isDynamicTimeFailureSuggested(this)){
		accepted = false;
		Tim->last_dt_accepted = false;
	}

	return Error;
}

#if 0
/**************************************************************************
   Copy the solution of the previous time interval to a vector
   temporarily
**************************************************************************/
void CRFProcessTH::StoreLastSolution(const int ty)
{
	long shift = 0;
	for (int i = 0; i < pcs_number_of_primary_nvals; i++)
	{
		const long number_of_nodes = num_nodes_p_var[i];
		for(long j = 0; j < number_of_nodes; j++)
			ARRAY[shift + j] = GetNodeValue(j, p_var_index[i] - ty);
		shift += number_of_nodes;
	}
}

void CRFProcessTH::ResetCouplingStep()
{
	long shift = 0;
	for (int i = 0; i < pcs_number_of_primary_nvals; i++)
	{
		const long number_of_nodes = num_nodes_p_var[i];
		for(long j = 0; j < number_of_nodes; j++)
			SetNodeValue(j, p_var_index[i], ARRAY[shift + j]);
		shift += number_of_nodes;
	}
}

void CRFProcessTH::InitializeNewtonSteps()
{
	for (int i = 0; i < pcs_number_of_primary_nvals; i++)
	{
		const int Col = p_var_index[i] - 1;
		const long number_of_nodes = num_nodes_p_var[i];
		for (int j = 0; j < number_of_nodes; j++)
			SetNodeValue(j, Col, 0.0);
	}
}
#endif

/**************************************************************************
   Update solution in Newton-Raphson procedure
**************************************************************************/
void CRFProcessTH::UpdateIterativeStep(const double damp)
{
	long shift = 0;
#if defined(NEW_EQS)
	const double* eqs_x = eqs_new->getX();
#elif !defined (USE_PETSC)
	const double* eqs_x = eqs->x;
#endif

	// x^k1 = x^k + dx
#if 0
	for (int i = 0; i < pcs_number_of_primary_nvals; i++)
	{
		const long number_of_nodes = num_nodes_p_var[i];
		const int ColIndex = p_var_index[i];
//		std::cout << std::endl;
//		for(long j = 0; j < number_of_nodes; j++) {
//			std::cout << eqs_x[j + shift] << " ";
//		}
		for (long j = 0; j < number_of_nodes; j++)
		{
#if defined(USE_PETSC)
			SetNodeValue(j, ColIndex, GetNodeValue(j,ColIndex) +
				damp * eqs_x[m_msh->Eqs2Global_NodeIndex[j]*pcs_number_of_primary_nvals + i]);
#else
			SetNodeValue(j, ColIndex, GetNodeValue(j, ColIndex) +  eqs_x[j + shift] * damp);
//			std::cout << GetNodeValue(j, ColIndex) << " ";
#endif
		}
//		std::cout << std::endl;
		shift += number_of_nodes;
	}
#endif
	const double p_fac = useMPa ? 1e6 : 1.0;
	const long number_of_nodes = num_nodes_p_var[0];
	const int ColIndex_p = p_var_index[0];
	if (m_num->petsc_split_fields) {
#if defined(USE_PETSC)
		Vec sub_x;
		std::vector<double> array_x(m_msh->getNumNodesGlobal());
		VecGetSubVector(eqs_new->x, eqs_new->vec_isg[0], &sub_x);
		eqs_new->getGlobalVectorArray(sub_x, &array_x[0]);
//		VecGetArray(sub_x, &array_x);
		//TODO collect ghost node values from other rank
		for (long j = 0; j < number_of_nodes; j++)
		{
			int k = m_msh->Eqs2Global_NodeIndex[j];
			SetNodeValue(j, ColIndex_p, GetNodeValue(j,ColIndex_p) + array_x[k] * damp * p_fac);
//			SetNodeValue(j, ColIndex_p, GetNodeValue(j,ColIndex_p) + array_x[j] * damp * p_fac);
		}
		VecRestoreSubVector(eqs_new->x, eqs_new->vec_isg[0], &sub_x);
		const int ColIndex_T = p_var_index[1];
		VecGetSubVector(eqs_new->x, eqs_new->vec_isg[1], &sub_x);
		eqs_new->getGlobalVectorArray(sub_x, &array_x[0]);
//		VecGetArray(sub_x, &array_x);
		for (long j = 0; j < number_of_nodes; j++)
		{
			int k = m_msh->Eqs2Global_NodeIndex[j];
			SetNodeValue(j, ColIndex_T, GetNodeValue(j,ColIndex_T) + array_x[k] * damp);
		}
		VecRestoreSubVector(eqs_new->x, eqs_new->vec_isg[1], &sub_x);

#endif
	} else {
		for (long j = 0; j < number_of_nodes; j++)
		{
#if defined(USE_PETSC)
			double dp =  eqs_x[m_msh->Eqs2Global_NodeIndex[j]*pcs_number_of_primary_nvals];
#else
			double dp = eqs_x[j] * damp;
	//			std::cout << GetNodeValue(j, ColIndex) << " ";
#endif
			SetNodeValue(j, ColIndex_p, GetNodeValue(j,ColIndex_p) + dp * damp * p_fac);
		}
		const int ColIndex_T = p_var_index[1];
		shift = number_of_nodes;
		for (long j = 0; j < number_of_nodes; j++)
		{
#if defined(USE_PETSC)
			double dT =  eqs_x[m_msh->Eqs2Global_NodeIndex[j]*pcs_number_of_primary_nvals + 1];
#else
			double dT = eqs_x[j+shift] * damp;
	//			std::cout << GetNodeValue(j, ColIndex) << " ";
#endif
			SetNodeValue(j, ColIndex_T, GetNodeValue(j,ColIndex_T) + dT * damp);
		}
	}

}

#if 0
/**************************************************************************
   Retrieve the solution from the temporary array
**************************************************************************/
void CRFProcessTH::RecoverSolution(const int ty)
{
	int i, j, idx;
	long number_of_nodes;
	int Colshift = 1;
	long shift = 0;
	double tem = 0.0;

	int start, end;

	start = 0;
	end = pcs_number_of_primary_nvals;

	for (i = start; i < end; i++)
	{
		number_of_nodes = num_nodes_p_var[i];
		idx =  p_var_index[i] - Colshift;
		for(j = 0; j < number_of_nodes; j++)
		{
			if(ty < 2)
			{
				if(ty == 1)
					tem =  GetNodeValue(j, idx);
				SetNodeValue(j, idx, ARRAY[shift + j]);
				if(ty == 1)
					ARRAY[shift + j] = tem;
			}
			else if(ty == 2)
			{
				tem = ARRAY[shift + j];
				ARRAY[shift + j] = GetNodeValue(j, idx);
				SetNodeValue(j,  idx,  tem);
			}
		}
		shift += number_of_nodes;
	}
}
#endif

