//  $Header$
//   Author: Toby Burnett
//  Project: geometry
//

#ifndef GEOMETRY_BOX_H
#define GEOMETRY_BOX_H

#include "geometry/Hexahedron.h"

/**
  Box is the same as the GEANT volume of the same name
  (except that full, rather than half, dimensions are specified)
  It is defined by six planes all perpendicular to axes
*/
class Box: public Hexahedron
{

 public:
   /// constructor: length along x, width along y, height along z
   Box( double length, double width, double height );

   /// resize Box: length along x, width along y, height along z
   void resize( double length, double width, double height );

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

