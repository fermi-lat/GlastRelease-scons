// ISpectrumFactory.h: interface for the ISpectrumFactory class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_ISPECTRUMFACTORY_H__10BF6F57_E416_4B69_A7CF_220E55281676__INCLUDED_)
#define AFX_ISPECTRUMFACTORY_H__10BF6F57_E416_4B69_A7CF_220E55281676__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000
#include "FluxSvc/Spectrum.h"

//!   This is an abstract base class for the SpectrumFactory template classes.

class ISpectrumFactory  
{
public:

	/// the only thing this does: make an associated Spectrum object
	virtual Spectrum* instantiate(const std::string& params)const=0;

	/// a dummy to force creation of the spectrum object.
	virtual void addRef()const=0;
};

#endif // !defined(AFX_ISPECTRUMFACTORY_H__10BF6F57_E416_4B69_A7CF_220E55281676__INCLUDED_)
