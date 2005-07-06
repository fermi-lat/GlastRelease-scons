
#ifndef RootConvert_Toolkit_H
#define RootConvert_Toolkit_H 1

#include <idents/VolumeIdentifier.h>
#include <commonRootData/idents/VolumeIdentifier.h>
#include <geometry/Point.h>
#include <geometry/Vector.h>
#include <CLHEP/Geometry/Point3D.h>
#include "TVector.h"
#include "TLorentzVector.h"
#include <string>

namespace RootPersistence {

/// low level conversions
idents::VolumeIdentifier convert( const VolumeIdentifier & ) ;
VolumeIdentifier convert( const idents::VolumeIdentifier & ) ;
HepLorentzVector convert( const TLorentzVector & ) ;
TLorentzVector convert( const HepLorentzVector & ) ;
HepPoint3D convert( const TVector3 & ) ;
TVector3 convert( const HepPoint3D & ) ;
Point convertPoint( const TVector3 & ) ;
TVector3 convertPoint( const Point & ) ;
Vector convertVector( const TVector3 & ) ;
TVector3 convertVector( const Vector & ) ;

/// low level comparisons
bool compare( Double_t, Double_t, const std::string & name ="unknown real" ) ;
bool compare( Int_t, Int_t, const std::string & name ="unknown integer" ) ;
bool compare( UInt_t, UInt_t, const std::string & name ="unknown integer" ) ;
bool compare( Bool_t, Bool_t, const std::string & name ="unknown bool" ) ;
bool compare( const std::string &, const std::string &, const std::string & name ="unknown vector3" ) ;
bool compare( const TVector3 &, const TVector3 &, const std::string & name ="unknown vector3" ) ;
bool compare( const TLorentzVector &, const TLorentzVector &, const std::string & name ="unknown lorentz vector" ) ;
bool compare( void *, void *, const std::string & name ="unknown pointer" ) ;

}

#endif

