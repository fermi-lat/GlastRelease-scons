/** @file MultipleScatteringFactory.h

   $Header$
   


*/
#ifndef MultipleScatteringFactory_h
#define  MultipleScatteringFactory_h

// from the G4 include monster
#include "G4VContinuousDiscreteProcess.hh"

namespace GlastMS {
    /** @class MultipleScatteringFactory
        @brief simple class to create a G4 process for Mscat.

    */
    class MultipleScatteringFactory {
    public:
        enum MScatType{ NATIVE, OLD32, LONGO};

        /** contstructor chooses which one
        @param type[NATIVE] get it from the installed G4 
        */
        MultipleScatteringFactory(MScatType type=NATIVE):m_type(type){}

        /**run the factory: this is a function operator
        */
     G4VContinuousDiscreteProcess* operator()();

    private:
       MScatType m_type;
    };

}// namespace

#endif
