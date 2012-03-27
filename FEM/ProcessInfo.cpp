/*
 * \file ProcessInfo.cpp
 *
 *  Created on: Sep 2, 2010
 *      Author: TF
 */

#include "rf_pcs.h"
#include <ProcessInfo.h>

ProcessInfo::ProcessInfo() :
	_pcs_type (FiniteElement::INVALID_PROCESS), _pcs_pv (FiniteElement::INVALID_PV), _pcs (NULL)
{}

ProcessInfo::ProcessInfo (FiniteElement::ProcessType pcs_type, FiniteElement::PrimaryVariable pcs_pv, CRFProcess* pcs,
						  FiniteElement::ErrorMethod nls_method, FiniteElement::ErrorMethod cpl_method) :
	_pcs_type (pcs_type), _pcs_pv (pcs_pv), _pcs (pcs), _pcs_nls_error_method(nls_method), _pcs_cpl_error_method(cpl_method)
{}

void ProcessInfo::setProcessType (FiniteElement::ProcessType pcs_type)
{
	_pcs_type = pcs_type;
}

void ProcessInfo::setProcessPrimaryVariable (FiniteElement::PrimaryVariable pcs_pv)
{
	_pcs_pv = pcs_pv;
}

void ProcessInfo::setProcess (CRFProcess* pcs)
{
	_pcs = pcs;
}

void ProcessInfo::setNonLinearErrorMethod (FiniteElement::ErrorMethod nls_method)
{
	_pcs_nls_error_method = nls_method;
}

void ProcessInfo::setCouplingErrorMethod (FiniteElement::ErrorMethod cpl_method)
{
	_pcs_cpl_error_method = cpl_method;
}

FiniteElement::ProcessType ProcessInfo::getProcessType () const
{
	return _pcs_type;
}

FiniteElement::PrimaryVariable ProcessInfo::getProcessPrimaryVariable () const
{
	return _pcs_pv;
}

CRFProcess* ProcessInfo::getProcess () const
{
	return _pcs;
}

FiniteElement::ErrorMethod ProcessInfo::getNonLinearErrorMethod () const
{
	return _pcs_nls_error_method;
}

FiniteElement::ErrorMethod ProcessInfo::getCouplingErrorMethod () const
{
	return _pcs_cpl_error_method;
}

ProcessInfo::~ProcessInfo()
{}
