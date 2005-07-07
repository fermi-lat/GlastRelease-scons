
#include "commonRootData/RootDataUtil.h"
#include "TVector.h"
#include "TLorentzVector.h"
#include "TString.h"
#include <iostream>

namespace rootdatautil {
    
void CheckBounds(long nIndex, long nLower, long nUpper) {
	if ((nIndex < nLower) || (nIndex > nUpper))
		throw(IndexOutOfBoundsException(nIndex, nLower, nUpper));
}

void CheckBufSize(unsigned long nBufSize, unsigned long nCorrectSize) {
	if (nBufSize != nCorrectSize)
		throw(BadBufferException(nBufSize, nCorrectSize));
}


//============================================================
// For tests, compare with an acceptable error for floats.
//============================================================

template <class T>
bool DirectCompareInRangeVals( T v1, T v2, const std::string & name ) {
  if (v1==v2) return true ;
  std::cout
    <<"Comparison ERROR for "<<name<<": "
    <<v1<<" != "<<v2
    <<std::endl ;
  return false ;
}

template <class T>
bool DirectCompareInRangeRefs( const T & v1, const T & v2, const std::string & name ) {
  if (v1==v2) return true ;
  std::cout
    <<"Comparison ERROR for "<<name<<": "
    <<v1<<" != "<<v2
    <<std::endl ;
  return false ;
}

bool CompareInRange( Double_t v1, Double_t v2, const std::string & name ) {
  const Double_t fudge=0.00001 ;
  if ((v1>=(v2-fudge))&&(v1<=(v2+fudge)) ){
    return true ;
  }
  std::cout
    <<"Comparison ERROR for "<<name<<": "
    <<v1<<" != "<<v2
    <<std::endl ;
  return false ;
}

bool CompareInRange( Int_t v1, Int_t v2, const std::string & name ) {
    return DirectCompareInRangeVals(v1,v2,name) ;
}
 
bool CompareInRange( UInt_t v1, UInt_t v2, const std::string & name ) {
    return DirectCompareInRangeVals(v1,v2,name) ;
}
 
bool CompareInRange( Bool_t v1, Bool_t v2, const std::string & name ) {
    return DirectCompareInRangeVals(v1,v2,name) ;
}
 
bool CompareInRange( const char * v1, const char * v2, const std::string & name ) {
    return DirectCompareInRangeRefs(std::string(v1),std::string(v2),name) ;
}
 
bool CompareInRange( const TString & v1, const TString & v2, const std::string & name ) {
    return DirectCompareInRangeRefs(std::string(v1.Data()),std::string(v2.Data()),name) ;
}
 
bool CompareInRange( const std::string & v1, const std::string & v2, const std::string & name ) {
    return DirectCompareInRangeRefs(v1,v2,name) ;
}
 
bool CompareInRange( const TVector3 & v1, const TVector3 & v2, const std::string & name ) {
    return CompareInRange(v1.X(),v2.X(),name+" X")
        && CompareInRange(v1.Y(),v2.Y(),name+" Y")
        && CompareInRange(v1.Z(),v2.Z(),name+" Z") ;
}
 
bool CompareInRange( const TLorentzVector & v1, const TLorentzVector & v2, const std::string & name ) {
    return CompareInRange(v1.X(),v2.X(),name+" X")
        && CompareInRange(v1.Y(),v2.Y(),name+" Y")
        && CompareInRange(v1.Z(),v2.Z(),name+" Z")
        && CompareInRange(v1.T(),v2.T(),name+" T") ;
}
 
bool CompareInRange( void * v1, void * v2, const std::string & name ) {
    return DirectCompareInRangeVals(v1,v2,name) ;
}
 
}
