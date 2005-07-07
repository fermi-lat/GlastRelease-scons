
#ifndef ROOTDATAUTIL_H
#define ROOTDATAUTIL_H

#include "Exception.h"

#include "Rtypes.h"
class TVector3 ;
class TLorentzVector ;

#include <string>

namespace rootdatautil {
    
    extern void CheckBounds(long nIndex, long nLower, long nUpper);
    extern void CheckBufSize(unsigned long nBufSize, unsigned long nCorrectSize); 

    // low level comparisons
    bool Compare( Double_t, Double_t, const std::string & name ="unknown real" ) ;
    bool Compare( Int_t, Int_t, const std::string & name ="unknown integer" ) ;
    bool Compare( UInt_t, UInt_t, const std::string & name ="unknown integer" ) ;
    bool Compare( Bool_t, Bool_t, const std::string & name ="unknown bool" ) ;
    bool Compare( const std::string &, const std::string &, const std::string & name ="unknown vector3" ) ;
    bool Compare( const TVector3 &, const TVector3 &, const std::string & name ="unknown vector3" ) ;
    bool Compare( const TLorentzVector &, const TLorentzVector &, const std::string & name ="unknown lorentz vector" ) ;
    bool Compare( void *, void *, const std::string & name ="unknown pointer" ) ;

};

#endif
