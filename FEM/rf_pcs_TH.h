#ifndef CRFPROCESSTH_H_
#define CRFPROCESSTH_H_

#include <vector>

#include "rf_pcs.h"


class CRFProcessTH : public CRFProcess
{
public:
	CRFProcessTH();
	virtual ~CRFProcessTH();

	void Initialization();
	virtual double Execute(int loop_process_number);
#if defined(USE_PETSC)
	virtual void setSolver( petsc_group::PETScLinearSolver *petsc_solver );
#endif

protected:
	void StoreLastSolution(const int ty=0);
	void ResetCouplingStep();
	void InitializeNewtonSteps();
	void RecoverSolution(const int ty = 0);
	void UpdateIterativeStep(const double damp);
	bool checkNRConvergence();

private:
	double error_k0;
	std::vector<double> ARRAY;

};

#endif
