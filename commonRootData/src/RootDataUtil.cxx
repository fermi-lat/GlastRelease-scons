
#include "commonRootData/RootDataUtil.h"
#include "TVector.h"
#include "TLorentzVector.h"
#include "TString.h"
#include "Riostream.h"

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

template <class T>
bool DelegatedCompareInRangeRefs( const T & v1, const T & v2, const std::string & name ) {
  return v1.CompareInRange(v2,name) ;
}

template <class T>
bool RealCompareInRangeVals( T v1, T v2, const std::string & name ) {
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

template <class StdCollection>
bool StdCollectionCompareInRange(
    const StdCollection & col1,
    const StdCollection & col2,
    const std::string & name ) {
    
    bool res = CompareInRange(col1.size(),col2.size(),name+" size") ;
    if (!res) return res ;
     
    typename StdCollection::const_iterator itr1 ;
    typename StdCollection::const_iterator itr2 ;
    for ( itr1 = col1.begin(), itr2 =col2.begin() ;
          ( itr1 != col1.end() ) && ( itr2 != col2.end() ) ;
          ++itr1, ++itr2 )
        res = res && CompareInRange(*itr1,*itr2,name+" element") ;
    
    return res ;
}


//============================================================
// Interface Implementation
//============================================================

bool CompareInRange( Float_t v1, Float_t v2, const std::string & name ) {
    return RealCompareInRangeVals(v1,v2,name) ;
}

bool CompareInRange( Double_t v1, Double_t v2, const std::string & name ) {
    return RealCompareInRangeVals(v1,v2,name) ;
}

bool CompareInRange( Short_t v1, Short_t v2, const std::string & name ) {
    return DirectCompareInRangeVals(v1,v2,name) ;
}

bool CompareInRange( Int_t v1, Int_t v2, const std::string & name ) {
    return DirectCompareInRangeVals(v1,v2,name) ;
}
 
bool CompareInRange( UInt_t v1, UInt_t v2, const std::string & name ) {
    return DirectCompareInRangeVals(v1,v2,name) ;
}

bool CompareInRange( ULong_t v1, ULong_t v2, const std::string & name ) {
    return DirectCompareInRangeVals(v1,v2,name) ;
}
 
bool CompareInRange( Bool_t v1, Bool_t v2, const std::string & name ) {
    return DirectCompareInRangeVals(v1,v2,name) ;
}
 
bool CompareInRange( Char_t v1, Char_t v2, const std::string & name ) {
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
 
bool CompareInRange( const TParameter<Double_t> & v1, const TParameter<Double_t> & v2, const std::string & name ) {
    return ( CompareInRange(v1.GetName(),v2.GetName(),name+" Parameter Name") &&
             CompareInRange(v1.GetVal(),v2.GetVal(),name+" Parameter Value") ) ;
}

bool CompareInRange( const AcdId & v1, const AcdId & v2, const std::string & name ) {
    return DelegatedCompareInRangeRefs(v1,v2,name) ;
}
 
bool CompareInRange( const CalXtalId & v1, const CalXtalId & v2, const std::string & name ) {
    return DelegatedCompareInRangeRefs(v1,v2,name) ;
}
 
bool CompareInRange( const CalXtalId::CalTrigMode & v1, const CalXtalId::CalTrigMode & v2, const std::string & name ) {
    return DirectCompareInRangeVals(v1,v2,name) ;
}

bool CompareInRange( const commonRootData::TkrId & v1, const commonRootData::TkrId & v2, const std::string & name ) {
    return DelegatedCompareInRangeRefs(v1,v2,name) ;
}

bool CompareInRange( const TowerId & v1, const TowerId & v2, const std::string & name ) {
    return DelegatedCompareInRangeRefs(v1,v2,name) ;
}

bool CompareInRange( const VolumeIdentifier & v1, const VolumeIdentifier & v2, const std::string & name ) {
    return DelegatedCompareInRangeRefs(v1,v2,name) ;
}

bool CompareInRange(
    const std::vector<Double_t> & col1,
    const std::vector<Double_t> & col2,
    const std::string & name ) {
    return StdCollectionCompareInRange(col1,col2,name) ;
}
bool CompareInRange(
    const std::vector<UInt_t> & col1,
    const std::vector<UInt_t> & col2,
    const std::string & name ) {
    return StdCollectionCompareInRange(col1,col2,name) ;
}

bool CompareInRange(
    const std::vector<AcdId> & col1,
    const std::vector<AcdId> & col2,
    const std::string & name ) {
    return StdCollectionCompareInRange(col1,col2,name) ;
}

}
