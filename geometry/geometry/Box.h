//  $Header$
//   Author: Toby Burnett
//  Project: geometry
//
//  Box is the same as the GEANT volume of the same name
//  (except that full, rather than half, dimensions are specified)
//  It is defined by six planes all perpendicular to axes

#ifndef BOX_H
#define BOX_H

#include "geometry/Hexahedron.h"


class Box: public Hexahedron
{

 public:
   Box( double length, double width, double height );
   // constructor: length along x, width along y, height along z

   void resize( double length, double width, double height );
   // resize Box: length along x, width along y, height along z

   double length()const;
   double width() const;
   double height()const;
   // return dimensions

   const Vector& xhat()const;
   const Vector& yhat()const;
   const Vector& zhat()const;
   // return references to the (original) x,y,z directions

 private:
   void constructSurfaces( double length, double width, double height );
   // make the 6 planes for this Box
};

#endif //BOX_H

