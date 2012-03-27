/*
 * \file ProcessInfo.h
 *
 *  Created on: Sep 2, 2010
 *      Author: TF
 */

#ifndef PROCESSINFO_H_
#define PROCESSINFO_H_

// FEM
#include "FEMEnums.h"
class CRFProcess;

/**
 * \brief Class ProcessInfo stores the process type,
 * an value for the primary variable of the process and
 * a pointer to the process object.
 */
class ProcessInfo
{
public:
	/**
	 * Default constructor, initializes pcs_type with ProcessType::INVALID_PROCESS,
	 * pcs_pv with PrimaryVariable::INVALID_PV and
	 * the pointer to the process with NULL. The user should set the values
	 * with the appropriate set methods.
	 * @return
	 */
	ProcessInfo();

	/**
	 * constructor initializing all attributes of the object with the given values
	 * @param pcs_type process type (\sa enum ProcessType)
	 * @param pcs_pv type of primary variable (\sa enum PrimaryVariable)
	 * @param pcs a pointer to the process
	 * @return
	 */
	ProcessInfo (FiniteElement::ProcessType pcs_type, FiniteElement::PrimaryVariable pcs_pv, CRFProcess* pcs,
		         FiniteElement::ErrorMethod nls_method, FiniteElement::ErrorMethod cpl_method);

	/**
	 * Sets the process type.
	 * @param pcs_type the process type, for valid values see enum ProcessType
	 */
	void setProcessType (FiniteElement::ProcessType pcs_type);

	/**
	 * Sets the value for the primary variable
	 * @param pcs_pv value for primary variable, possible values are documented in enum PrimaryVariable
	 */
	void setProcessPrimaryVariable (FiniteElement::PrimaryVariable pcs_pv);

	/**
	 * Sets the value for the pointer to an object of class CRFProcess.
	 * @param pcs the pointer to an object of class CRFProcess
	 */
	void setProcess (CRFProcess* pcs);

	/**
	 * Sets the process non-linear error method
	 */
	void setNonLinearErrorMethod (FiniteElement::ErrorMethod nls_method);

	/**
	 * Sets the process coupling error method
	 */
	void setCouplingErrorMethod (FiniteElement::ErrorMethod cpl_method);

	/**
	 * Get the process type.
	 * @return the process type
	 */
	FiniteElement::ProcessType getProcessType () const;

	/**
	 * Get the primary variable of the process.
	 * @return the primary variable of the process
	 */
	FiniteElement::PrimaryVariable getProcessPrimaryVariable () const;

	/**
	 * Get a pointer to an object of type CRFProcess.
	 * @return a pointer to an object of type CRFProcess
	 */
	CRFProcess* getProcess () const;

	/**
	 * Gets the process non-linear error method
	 */
	FiniteElement::ErrorMethod getNonLinearErrorMethod () const;

	/**
	 * Gets the process coupling error method
	 */
	FiniteElement::ErrorMethod getCouplingErrorMethod () const;

	virtual ~ProcessInfo();

protected:
	/**
	 * process type, see enum ProcessType for valid values
	 */
	FiniteElement::ProcessType _pcs_type;
	/**
	 * the primary variable of the process, see enum PrimaryVariable for valid values
	 */
	FiniteElement::PrimaryVariable _pcs_pv;
	/**
	 * the non-linear error method of the process
	 */
	FiniteElement::ErrorMethod _pcs_nls_error_method;
	/**
	 * the coupling error method of the process
	 */
	FiniteElement::ErrorMethod _pcs_cpl_error_method;

	/**
	 * pointer to the object of class CRFProcess
	 */
	CRFProcess* _pcs;
};
#endif                                            /* PROCESSINFO_H_ */
