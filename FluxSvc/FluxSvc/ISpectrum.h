// ISpectrum.h: interface for the ISpectrum class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_ISPECTRUM_H__9D00078B_7121_473C_813A_827B2F4126A5__INCLUDED_)
#define AFX_ISPECTRUM_H__9D00078B_7121_473C_813A_827B2F4126A5__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <string>
#include <utility>
#include <vector>
// CLHEP
#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/RandGeneral.h"
#include "CLHEP/Random/RandExponential.h"
#include "CLHEP/Random/RandFlat.h"

//forward declaration
class HepRandomEngine;


class ISpectrum  
{


//!  Class for holding function definitions of Spectrums...
//!  an abstract base class
public:
    
    
    virtual float operator()(float /*r*/)const=0;
    
    /// subclasses need to specify correct particle type
    virtual const char * particleName()const=0;
    
    /// calculate the flux, particles/m^2/sr. (default zero)
    virtual double    flux (double time ) const=0;
    
    
    /// calcualte effective solid angle  (default zero)
    virtual double solidAngle()const{return 6.;}
    
    /// return a title describing the spectrum	
    virtual std::string title()const=0;
    
    /// a randomized interval to the next event - default is 1/rate()
    virtual double interval (double time)=0;

    virtual std::pair<float,float> dir(float energy)const = 0;
    
    /// new interface for Hirosima classes
    virtual double energySrc(HepRandomEngine* engine, double time=0)=0;
    virtual std::pair<double,double> dir(double energy, HepRandomEngine* engine)=0;

};


#endif // !defined(AFX_ISPECTRUM_H__9D00078B_7121_473C_813A_827B2F4126A5__INCLUDED_)
