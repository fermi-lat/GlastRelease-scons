// $Id$
//
//

#ifndef GEOMETRY_GEOMOBJECT_H
#define GEOMETRY_GEOMOBJECT_H


#include <iostream>

extern void WARNING(const char *);
extern void FATAL(const char *);
// handle error messages: should make it create exceptions?

class CoordTransform;
class Vector;
class Command;

/// Pure base class for all Geometric objects. Establishes the ability to transform.
class GeomObject 
{
 public:
   virtual ~GeomObject(){}

   // Implicit copy constructor
   // Implicit assignment operator

   /// virtual method to perform transformations that must be implemented by subclasses
   virtual GeomObject& transform(const CoordTransform&)=0;

   /// gives subclasses ability to just move
   GeomObject& move(const Vector&);
   GeomObject& moveX(double);
   GeomObject& moveY(double);
   GeomObject& moveZ(double);

   /// Rotation around the given axis: return pointer to self for convenience
   GeomObject & rotateX(double);
   GeomObject & rotateY(double);
   GeomObject & rotateZ(double);

    /// Rotation around a specified vector.
    GeomObject & rotate(double angle, const Vector & axis);

    /// makes name available to subclasses using the type info facility
    virtual const char*  nameOf()const;

    /// print info about this object
    virtual void printOn(std::ostream& out=std::cout)const;
 


protected:
    GeomObject(){}

};

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~inlines~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

inline std::ostream& operator<<(std::ostream& os, const GeomObject& g){
   g.printOn(os); return os;}

#endif

