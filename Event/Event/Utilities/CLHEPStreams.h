// $Header$
#ifndef LHCBEVENT_CLHEPSTREAMS_H
#define LHCBEVENT_CLHEPSTREAMS_H 1


// Include files
#include "GaudiKernel/StreamBuffer.h"
#include "Event/TopLevel/Definitions.h"


//------------------------------------------------------------------------------
//
// Description: Streams operators of CLHEP classes used in Event
//              (used in serialize() methods)
//
// CLHEPStreams.h defines additional oprators to used classes from
// <A HREF="http://wwwinfo.cern.ch/asd/lhc++/clhep/manual/RefGuide/index.html">CLHEP</A>
//              
// Author:      Pavel Binko
//
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
// <A HREF="http://wwwinfo.cern.ch/asd/lhc++/clhep/manual/RefGuide/Geometry/HepPoint3D.html">class HepPoint3D</A>
#include "CLHEP/Geometry/Point3D.h"
#include "CLHEP/config/CLHEP.h"
// Hack for CLHEP 1.9.2.2
#ifndef HepPoint3D
typedef HepGeom::Point3D<double> HepPoint3D;
#endif
//------------------------------------------------------------------------------
// Output operator
inline StreamBuffer& operator<< ( StreamBuffer& s, const HepPoint3D& value)    {
  return s  << value.x() << value.y() << value.z();
}
// Input operator
inline StreamBuffer& operator>> ( StreamBuffer& s, HepPoint3D& value )         {
  HepDouble   x, y, z;
  s >> x >> y >> z;
  value.setX(x);
  value.setY(y);
  value.setZ(z);
  return s;
}


//------------------------------------------------------------------------------
// <A HREF="http://wwwinfo.cern.ch/asd/lhc++/clhep/manual/RefGuide/Geometry/HepVector3D.html">class HepVector3D</A>
#include "CLHEP/Geometry/Vector3D.h"
// Hack for CLHEP 1.9.2.2
#ifndef HepVector3D
typedef HepGeom::Vector3D<double> HepVector3D;
#endif
//------------------------------------------------------------------------------
// Output operator
inline StreamBuffer& operator<< ( StreamBuffer& s, const HepVector3D& value)   {
  return s  << value.x() << value.y() << value.z();
}
// Input operator
inline StreamBuffer& operator>> ( StreamBuffer& s, HepVector3D& value )        {
  HepDouble   x, y, z;
  s  >> x >> y >> z;
  value.setX(x);
  value.setY(y);
  value.setZ(z);
  return s;
}


//------------------------------------------------------------------------------
// <A HREF="http://wwwinfo.cern.ch/asd/lhc++/clhep/manual/RefGuide/Vector/HepLorentzVector.html">class HepLorentzVector</A>
#include "CLHEP/Vector/LorentzVector.h"
//------------------------------------------------------------------------------
// Output operator
inline StreamBuffer& operator<< ( StreamBuffer& s, const CLHEP::HepLorentzVector& value) {
  return s  << HepDouble(value.px()) 
            << HepDouble(value.py())
            << HepDouble(value.pz())
            << HepDouble(value.e());
}
// Input operator
inline StreamBuffer& operator>> ( StreamBuffer& s, CLHEP::HepLorentzVector& value )   {
  HepDouble   px, py, pz, E;
  s  >> px >> py >> pz >> E;
  value.setPx(px);
  value.setPy(py);
  value.setPz(pz);
  value.setE(E);
  return s;
}


//------------------------------------------------------------------------------
// <A HREF="http://wwwinfo.cern.ch/asd/lhc++/clhep/manual/RefGuide/Matrix/HepSymMatrix.html">class SymMatrix</A>
#include "CLHEP/Matrix/SymMatrix.h"
//------------------------------------------------------------------------------
// Output operator
inline StreamBuffer& operator<< ( StreamBuffer& s, const CLHEP::HepSymMatrix& value ) {
  int   nrow = value.num_row();
  for( int i=1; i<=nrow; i++ ) {
    for( int j=1; j<=i; j++ ) {
      s << value(i,j);
    }
  }
  return s;
}
// Input operator
inline StreamBuffer& operator>> ( StreamBuffer& s, CLHEP::HepSymMatrix& value )       {
  int   nrow = value.num_row();
  for( int i=1; i<=nrow; i++ ) {
    for( int j=1; j<=i; j++ ) {
      s >> value(i,j);
    }
  }
  return s;
}


#endif    // LHCBEVENT_CLHEPSTREAMS_H
