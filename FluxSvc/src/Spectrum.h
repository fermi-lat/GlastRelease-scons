// $Header$
//
//
// Spectrum: base class for energy spectrum objects
// SimpleSpectrum: define a particle and spectral index
//

#ifndef GLAST_SPECTRUM_H
#define GLAST_SPECTRUM_H

#include <string>
#include <utility>
#include <vector>
#include "FluxSvc/ISpectrum.h"
//forward declaration
class HepRandomEngine;


//!  Class for holding function definitions of Spectrums - i.e. HeSpectrum, SimpleSpectrum, etc...
//!  Basically an abstract base class for these classes.
class Spectrum : public ISpectrum {
public:
    
    //    class  Direction : public std::pair<float,float> {
    //    public:
    //        double costh()const{return this.first;}
    //        double phi()const {return this.second;}
    //    };
    
    virtual float operator()(float /*r*/)const{return 0;};
    // this is all that these objects do. Must be virtual for
    // polymorphism
    // returns kinetic energy for random number r in [0,1). 
    // requried that it be monatonic
    // NB. Default is to return zero, an indicator that the actual Spectrum object
    //     implements a method that makes direct use of the random generator
    
    //virtual double calculate_rate (double old_rate) = 0;
    
    /// subclasses need to specify correct particle type
    virtual const char * particleName()const=0;
    
    /// calculate the flux, particles/m^2/sr. (default zero)
    virtual double    flux (double time ) const;
    
    
    /// calcualte effective solid angle  (default zero)
    virtual double solidAngle()const;
    
    /// return a title describing the spectrum	
    virtual std::string title()const=0;
    
    /// inverse of the operator: given a KE, return number in [0,1)
    /// for choosing limits
    float fraction(float energy);
    
    virtual ~Spectrum();
    
    /// a randomized interval to the next event - default is 1/rate()
    virtual double interval (double time);
    
    virtual std::pair<float,float> dir(float energy)const;
    
    /// new interface for Hirosima classes
    virtual double energySrc(HepRandomEngine* engine, double time=0);
    virtual std::pair<double,double> dir(double energy, HepRandomEngine* engine);
    
    
    
    
protected:
    Spectrum(const std::vector<float>& /*params*/){};
    Spectrum(/* double lat = 0, double lon = 0, double time=0*/) 
        /*: m_lat(0), m_lon(0), m_time(0)*/{}
        // all constructors protected to ensure an abstract class
        
        virtual void parseParamList(std::string input, std::vector<float>& output) const;
    
    double    m_lat, m_lon;   // latitude and longitudinal coordinates
    double m_currentInterval; // so we only find the interval for each particle once.
    
};

#endif    
