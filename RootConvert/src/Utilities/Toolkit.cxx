
#include <RootConvert/Utilities/Toolkit.h>
#include <iostream>

namespace RootPersistence {

//============================================================
// idents conversions
//============================================================

idents::AcdId convert( const AcdId & rootId ) {
  //throw "idents::AcdId convert( const AcdId & ) NOT IMPLEMENTED" ;    
    return idents::AcdId(rootId.getId(2));
}


AcdId convert( const idents::AcdId & tdsId ) {
  //throw "AcdId convert( const idents::AcdId & ) NOT IMPLEMENTED" ;
    AcdId root;  root.setId( tdsId.id() );
    return root;
}


idents::CalXtalId convert( const CalXtalId & rootId ) {
    return idents::CalXtalId(rootId.getTower(),rootId.getLayer(),rootId.getColumn()) ;
}


CalXtalId convert( const idents::CalXtalId & tdsId ) {
    return CalXtalId(tdsId.getTower(),tdsId.getLayer(),tdsId.getColumn()) ;
}


idents::TkrId convert( const commonRootData::TkrId & ) {
    throw "idents::TkrId convert( const commonRootData::TkrId & ) NOT IMPLEMENTED" ;
    return idents::TkrId() ;
}


commonRootData::TkrId convert( const idents::TkrId & ) {
    throw "commonRootData::TkrId convert( const idents::TkrId & ) NOT IMPLEMENTED" ;
    return commonRootData::TkrId() ;
}


idents::TowerId convert( const TowerId & ) {
    throw "idents::TowerId convert( const TowerId & ) NOT IMPLEMENTED" ;
    return idents::TowerId() ;
}


TowerId convert( const idents::TowerId & ) {
    throw "TowerId convert( const idents::TowerId & ) NOT IMPLEMENTED" ;
    return TowerId() ;
}

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
    //           in the ROOT class.  The idents::VolumeIdentifier class stores the
    //           data in one 64 bit word.  We must convert from the two 32 bit words
    //           into the 64 bit word.  We perform the conversion by iterating over
    //           all of the ids in the ROOT VolumeIdentifier and appending them to
    //           the TDS idents::VolumeIdentifier.
    // Input:  ROOT VolumeIdentifier
    // Ouput:  idents::VolumeIdentifier
    
    idents::VolumeIdentifier tdsVolId ;
    int index;
    for (index = 0; index < rootVolId.size(); index++) {
        tdsVolId.append(rootVolId.operator [](index));
    }
    return tdsVolId ;
    
}


//============================================================
// Conversion math
//============================================================
 typedef HepGeom::Point3D<double> HepPoint3D;


HepPoint3D convert( const TVector3 & v ) {
    return HepPoint3D(v.X(),v.Y(),v.Z()) ;
}

TVector3 convert( const HepPoint3D & v ) {
    return TVector3(v.x(),v.y(),v.z()) ;
}

CLHEP::Hep3Vector convert3vector( const TVector3 & v ) {
    return CLHEP::Hep3Vector(v.X(),v.Y(),v.Z()) ;
}

TVector3 convert3vector( const CLHEP::Hep3Vector & v ) {
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

CLHEP::HepLorentzVector convert( const TLorentzVector & v ) {
    return CLHEP::HepLorentzVector(v.X(),v.Y(),v.Z(),v.T()) ;
}

TLorentzVector convert( const CLHEP::HepLorentzVector & v ) {
    return TLorentzVector(v.x(),v.y(),v.z(),v.t()) ;
}

}
