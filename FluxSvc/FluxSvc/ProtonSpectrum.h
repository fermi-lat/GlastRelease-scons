// $Id$
/* $History: ProtonSpectrum.h $
* 
* *****************  Version 2  *****************
* User: Burnett      Date: 29-03-98   Time: 15:54
* Updated in $/glastsim/flux
* Change Flux definition so that subclasses must implement a descriptive
* title
*/

//
// return the cosmic ray proton spectrum appropriate for low 
// inclination Earth orbit
//

#ifndef PROTON_SPECTRUM_H
#define PROTON_SPECTRUM_H

#include "FluxSvc/Spectrum.h"


//! This is the cosmic ray proton spectrum appropriate for low 
//! inclination Earth orbit.
class ProtonSpectrum : public Spectrum
{
public:
    ProtonSpectrum(){};
    
    virtual float operator()(float)const;
    virtual double calculate_rate(double old_rate);
    virtual const char* particleName()const;
    virtual std::string title()const;
    
private:
    static int ProtonSpectrum::prolocate(struct tbl *data, int n, double x);
    static double protonenergy(float r);
    
};

#endif  
