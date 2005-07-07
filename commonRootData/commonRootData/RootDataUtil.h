
#ifndef ROOTDATAUTIL_H
#define ROOTDATAUTIL_H

#include "Exception.h"

#include "Rtypes.h"
class TVector3 ;
class TLorentzVector ;
class TString ;

#include <string>

namespace rootdatautil {
    
    extern void CheckBounds(long nIndex, long nLower, long nUpper);
    extern void CheckBufSize(unsigned long nBufSize, unsigned long nCorrectSize); 

    // low level comparisons
    bool CompareInRange( Double_t, Double_t, const std::string & name ="unknown real" ) ;
    bool CompareInRange( Int_t, Int_t, const std::string & name ="unknown integer" ) ;
    bool CompareInRange( UInt_t, UInt_t, const std::string & name ="unknown integer" ) ;
    bool CompareInRange( Bool_t, Bool_t, const std::string & name ="unknown bool" ) ;
    bool CompareInRange( const char *, const char *, const std::string & name ="unknown vector3" ) ;
    bool CompareInRange( const TString &, const TString &, const std::string & name ="unknown vector3" ) ;
    bool CompareInRange( const std::string &, const std::string &, const std::string & name ="unknown vector3" ) ;
    bool CompareInRange( const TVector3 &, const TVector3 &, const std::string & name ="unknown vector3" ) ;
    bool CompareInRange( const TLorentzVector &, const TLorentzVector &, const std::string & name ="unknown lorentz vector" ) ;
    bool CompareInRange( void *, void *, const std::string & name ="unknown pointer" ) ;

};

#endif
