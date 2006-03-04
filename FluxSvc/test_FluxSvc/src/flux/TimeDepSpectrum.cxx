#include "TimeDepSpectrum.h"
// define a factory for anonomous instantiation
#include "../../SpectrumFactory.h"

static SpectrumFactory<TimeDepSpectrum> factory;
const ISpectrumFactory& TimeDepSpectrumFactory = factory;

 /// return total flux 
/*double TimeDepSpectrum::flux() const{
       double time = (GPS::instance()->time());
       if(time<=2.) return 1;
       return 0;
   }
  */ 
   ///return flux, given a time
   double TimeDepSpectrum::flux(double time)const{
       if(time<=0.04) {return 1;}
       return 0;
   }
