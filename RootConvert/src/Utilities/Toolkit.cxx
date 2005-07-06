
#include <RootConvert/Utilities/Toolkit.h>
#include <iostream>

namespace RootPersistence {

//============================================================
// Conversion
//============================================================

VolumeIdentifier convert( const idents::VolumeIdentifier & tdsVolId ) {
    
    // Purpose and Method:  We must store the volume ids as two 32 bit UInt_t
    //     in the ROOT class.  Hence, we must convert the 64 bit representation
    //     used in the idents::VolumeIdentifier class into two 32 bit UInt_t.
    //     To perform the conversion, we iterate over all the ids in the TDS
    //     version of the idents::VolumeIdentifier and append each to the ROOT
    //     VolumeIdentifier
    
    VolumeIdentifier rootVolId ;
    int index;
    rootVolId.Clear();
    for (index = 0; index < tdsVolId.size(); index++) {
        rootVolId.append(tdsVolId.operator [](index));
    }
    return rootVolId ;
}
 
idents::VolumeIdentifier convert( const VolumeIdentifier & rootVolId ) {

    // Purpose and Method:  We must store the volume ids as two 32 bit UInt_t
    //	   in the ROOT class.  The idents::VolumeIdentifier class stores the
    //	   data in one 64 bit word.  We must convert from the two 32 bit words
    //	   into the 64 bit word.  We perform the conversion by iterating over
    //	   all of the ids in the ROOT VolumeIdentifier and appending them to
    //	   the TDS idents::VolumeIdentifier.
    // Input:  ROOT VolumeIdentifier
    // Ouput:  idents::VolumeIdentifier
    
    idents::VolumeIdentifier tdsVolId ;
    int index;
    for (index = 0; index < rootVolId.size(); index++) {
        tdsVolId.append(rootVolId.operator [](index));
    }
    return tdsVolId ;
    
}

HepPoint3D convert( const TVector3 & v ) {
    return HepPoint3D(v.X(),v.Y(),v.Z()) ;
}

TVector3 convert( const HepPoint3D & v ) {
    return TVector3(v.x(),v.y(),v.z()) ;
}

Point convertPoint( const TVector3 & v ) {
    return Point(v.X(),v.Y(),v.Z()) ;
}

TVector3 convertPoint( const Point & p ) {
    return TVector3(p.getX(),p.getY(),p.getZ()) ;
}

Vector convertVector( const TVector3 & v ) {
    return Vector(v.X(),v.Y(),v.Z()) ;
}

TVector3 convertVector( const Vector & v ) {
    return TVector3(v.getX(),v.getY(),v.getZ()) ;
}

HepLorentzVector convert( const TLorentzVector & v ) {
    return HepLorentzVector(v.X(),v.Y(),v.Z(),v.T()) ;
}

TLorentzVector convert( const HepLorentzVector & v ) {
    return TLorentzVector(v.x(),v.y(),v.z(),v.t()) ;
}

//============================================================
// Comparison
//============================================================

template <class T>
bool direct_compare_vals( T v1, T v2, const std::string & name ) {
  if (v1==v2) return true ;
  std::cout
    <<"Comparison failed for "<<name<<": "
    <<v1<<" != "<<v2
    <<std::endl ;
  return false ;
}

template <class T>
bool direct_compare_refs( const T & v1, const T & v2, const std::string & name ) {
  if (v1==v2) return true ;
  std::cout
    <<"Comparison failed for "<<name<<": "
    <<v1<<" != "<<v2
    <<std::endl ;
  return false ;
}

bool compare( Double_t v1, Double_t v2, const std::string & name ) {
  const Double_t fudge=0.00001 ;
  if ((v1>=(v2-fudge))&&(v1<=(v2+fudge)) ){
    return true ;
  }
  std::cout
    <<"Comparison failed for "<<name<<": "
    <<v1<<" != "<<v2
    <<std::endl ;
  return false ;
}

bool compare( Int_t v1, Int_t v2, const std::string & name ) {
    return direct_compare_vals(v1,v2,name) ;
}
 
bool compare( UInt_t v1, UInt_t v2, const std::string & name ) {
    return direct_compare_vals(v1,v2,name) ;
}
 
bool compare( Bool_t v1, Bool_t v2, const std::string & name ) {
    return direct_compare_vals(v1,v2,name) ;
}
 
bool compare( const std::string & v1, const std::string & v2, const std::string & name ) {
    return direct_compare_refs(v1,v2,name) ;
}
 
bool compare( const TVector3 & v1, const TVector3 & v2, const std::string & name ) {
    return compare(v1.X(),v2.X(),name+" X")
        && compare(v1.Y(),v2.Y(),name+" Y")
        && compare(v1.Z(),v2.Z(),name+" Z") ;
}
 
bool compare( const TLorentzVector & v1, const TLorentzVector & v2, const std::string & name ) {
    return compare(v1.X(),v2.X(),name+" X")
        && compare(v1.Y(),v2.Y(),name+" Y")
        && compare(v1.Z(),v2.Z(),name+" Z")
        && compare(v1.T(),v2.T(),name+" T") ;
}
 
bool compare( void * v1, void * v2, const std::string & name ) {
    return direct_compare_vals(v1,v2,name) ;
}
 

}
