// $Header$
//
// Simulate the Earth albedo neutron and photon spectra.
// Choose origins uniformly in a plane, and isotropic angles
//
#ifndef ALBEDO_SPECTRUM_H
#define ALBEDO_SPECTRUM_H

#include "FluxSvc/Spectrum.h"

namespace xml {
  class Element;
}

//! Simulate the Earth albedo neutron and photon spectra.
//! Choose origins uniformly in a plane, and isotropic angles
class AlbedoSpectrum : public Spectrum
{
public:
    AlbedoSpectrum(){
		std::string str="gamma";
		AlbedoSpectrum("gamma");
	}

    /// require attribute name="gamma" or name="neutron" 
    AlbedoSpectrum(const xml::Element& xelem);

    /// constructor. must specify name ("photon" or "n")
    AlbedoSpectrum(const std::string& name);

    virtual double calculate_rate(double old_rate);


    /// return energy based on random
    virtual float operator()(float r)const;
    
    virtual const char* particleName()const;
    virtual std::string title()const;
    
    enum Type{photon, neutron};
    
    void setPosition ( float lat, float lon ){}
    
private:

    ///photon or neutron (or...?)
    Type type_; 
    
};

#endif     
