// SpectrumFactory.h: interface for the SpectrumFactory class.
//
//////////////////////////////////////////////////////////////////////
// $Header$

#if !defined(AFX_SPECTRUMFACTORY_H__211C2F25_9111_44B9_B357_0762789222AF__INCLUDED_)
#define AFX_SPECTRUMFACTORY_H__211C2F25_9111_44B9_B357_0762789222AF__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "FluxSvc/ISpectrumFactory.h"
#include "FluxSvc/SpectrumFactoryTable.h"
#include <typeinfo>
#include <vector>


//!  Template class designed to hold method by which to polymorphically instantiate new Spectrum objects

template <class T> class SpectrumFactory : public ISpectrumFactory 
{
public:
    
    SpectrumFactory(){
        //Get class name using RTTI:
        std::string classname = typeid(T).name();
        int s = classname.find_first_of("class");
        if( s==0 ) s=6; //found it
        else s =classname.find_first_not_of("0123456789");
        classname = classname.substr(s);
        SpectrumFactoryTable::instance()->addFactory(classname, this); 
    }
    //! return a new Spectrum object
    virtual Spectrum* instantiate(const std::string& params) const{return new T(params);}
    
    //! dummy to follow Gaudi model
    virtual void addRef()const{}
};

#endif // !defined(AFX_SPECTRUMFACTORY_H__211C2F25_9111_44B9_B357_0762789222AF__INCLUDED_)
