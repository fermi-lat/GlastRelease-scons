
#ifndef RootConvert_Toolkit_H
#define RootConvert_Toolkit_H 1

#include <idents/AcdId.h>
#include <idents/CalXtalId.h>
#include <idents/TkrId.h>
#include <idents/TowerId.h>
#include <idents/VolumeIdentifier.h>
#include <commonRootData/idents/AcdId.h>
#include <commonRootData/idents/CalXtalId.h>
#include <commonRootData/idents/TkrId.h>
#include <commonRootData/idents/TowerId.h>
#include <commonRootData/idents/VolumeIdentifier.h>

#include <geometry/Point.h>
#include <geometry/Vector.h>
#include <CLHEP/Geometry/Point3D.h>
#include <CLHEP/Vector/LorentzVector.h>
#include "TVector.h"
#include "TLorentzVector.h"
#include <string>

// Hack for CLHEP 1.9.2.2
#ifndef HepPoint3D
typedef HepGeom::Point3D<double> HepPoint3D;
#endif

namespace RootPersistence {

idents::AcdId convert( const AcdId & ) ;
AcdId convert( const idents::AcdId & ) ;
idents::CalXtalId convert( const CalXtalId & ) ;
CalXtalId convert( const idents::CalXtalId & ) ;
idents::TkrId convert( const commonRootData::TkrId & ) ;
commonRootData::TkrId convert( const idents::TkrId & ) ;
idents::TowerId convert( const TowerId & ) ;
TowerId convert( const idents::TowerId & ) ;
idents::VolumeIdentifier convert( const VolumeIdentifier & ) ;
VolumeIdentifier convert( const idents::VolumeIdentifier & ) ;
CLHEP::HepLorentzVector convert( const TLorentzVector & ) ;
TLorentzVector convert( const CLHEP::HepLorentzVector & ) ;
HepPoint3D convert( const TVector3 & ) ;
TVector3 convert( const HepPoint3D & ) ;
CLHEP::Hep3Vector convert3vector( const TVector3 & ) ;
TVector3 convert3vector( const CLHEP::Hep3Vector & ) ;
Point convertPoint( const TVector3 & ) ;
TVector3 convertPoint( const Point & ) ;
Vector convertVector( const TVector3 & ) ;
TVector3 convertVector( const Vector & ) ;

}

#endif

