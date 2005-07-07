
#include "commonRootData/RootDataUtil.h"
#include "TVector.h"
#include "TLorentzVector.h"
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
bool DirectCompareVals( T v1, T v2, const std::string & name ) {
  if (v1==v2) return true ;
  std::cout
    <<"Comparison ERROR for "<<name<<": "
    <<v1<<" != "<<v2
    <<std::endl ;
  return false ;
}

template <class T>
bool DirectCompareRefs( const T & v1, const T & v2, const std::string & name ) {
  if (v1==v2) return true ;
  std::cout
    <<"Comparison ERROR for "<<name<<": "
    <<v1<<" != "<<v2
    <<std::endl ;
  return false ;
}

bool Compare( Double_t v1, Double_t v2, const std::string & name ) {
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

bool Compare( Int_t v1, Int_t v2, const std::string & name ) {
    return DirectCompareVals(v1,v2,name) ;
}
 
bool Compare( UInt_t v1, UInt_t v2, const std::string & name ) {
    return DirectCompareVals(v1,v2,name) ;
}
 
bool Compare( Bool_t v1, Bool_t v2, const std::string & name ) {
    return DirectCompareVals(v1,v2,name) ;
}
 
bool Compare( const std::string & v1, const std::string & v2, const std::string & name ) {
    return DirectCompareRefs(v1,v2,name) ;
}
 
bool Compare( const TVector3 & v1, const TVector3 & v2, const std::string & name ) {
    return Compare(v1.X(),v2.X(),name+" X")
        && Compare(v1.Y(),v2.Y(),name+" Y")
        && Compare(v1.Z(),v2.Z(),name+" Z") ;
}
 
bool Compare( const TLorentzVector & v1, const TLorentzVector & v2, const std::string & name ) {
    return Compare(v1.X(),v2.X(),name+" X")
        && Compare(v1.Y(),v2.Y(),name+" Y")
        && Compare(v1.Z(),v2.Z(),name+" Z")
        && Compare(v1.T(),v2.T(),name+" T") ;
}
 
bool Compare( void * v1, void * v2, const std::string & name ) {
    return DirectCompareVals(v1,v2,name) ;
}
 
}