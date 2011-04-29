/** @file MultipleScatteringFactory.cxx

   $Header$
   


*/
#include "GlastMS/MultipleScatteringFactory.h"

//#include "G4eMultipleScattering.hh"
// To be changed (to have MCS per particle)

// from this package
#include "GlastMS/G4MultipleScattering.h"
#include "GlastMS/G4MultipleScatteringFL.h"

namespace GlastMS{

G4VContinuousDiscreteProcess* MultipleScatteringFactory::operator()(){
        if( m_type==NATIVE){
            // the G4 version currently out there
            //	return new ::G4eMultipleScattering();
            return new GlastMS::G4MultipleScattering();
        }
        else if (m_type == OLD32)
        {
            // the local copy
            return new GlastMS::G4MultipleScattering();
        }
        else
        {
            // the local copy
            return new GlastMS::G4MultipleScatteringFL();
        }
    }
}
