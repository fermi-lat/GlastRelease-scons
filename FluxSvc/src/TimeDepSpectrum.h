#ifndef TimeDepSpectrum_H
#define TimeDepSpectrum_H
// $Heading:$
//
#include "FluxSvc/Spectrum.h"
#include "facilities/Observer.h"
#include <vector>
#include <string>
#include <map>
#include "GPS.h"

//!  Spectrum that reads its differential spectrum from a table
class TimeDepSpectrum : public Spectrum
{
public:
   /// params is the filename to read
    TimeDepSpectrum(const std::string& params){}

   
   /// return total flux 
//   virtual double flux() const;
   
   ///return flux, given a time
   double flux(double time)const;
   
   /// sample a single particle energy from the spectrum
   inline float operator() (float)const{return (float)0.1;}
   
   
   inline std::string title() const {return "TimeDepSpectrum";}
   inline const char * particleName() const {return "e";}
   inline  const char * nameOf() const {return "TimeDepSpectrum";}
   //   use default destructor, copy constructor, and assignment op.
   
   
private:
   
};



#endif