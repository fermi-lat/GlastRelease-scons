// $Id$
//
//
// Pure base class for all Geometric objects. Establishes the ability
// to transform.
//
// Inherits from ArvePlottable to provide some default display capability

#ifndef GEOMOBJECT_H
#define GEOMOBJECT_H


#include <iostream>

extern void WARNING(const char *);
extern void FATAL(const char *);
// handle error messages: should make it create exceptions?

class CoordTransform;
class Vector;
class Command;
// just declare these, so the .cxx files need to include them

class GeomObject //: public ArvePlottable
{
 public:
   virtual ~GeomObject(){}

   // Implicit copy constructor
   // Implicit assignment operator

   virtual GeomObject& transform(const CoordTransform&)=0;
   // virtual method to perform transformations that must be
   // implemented by subclasses

   GeomObject& move(const Vector&);
   GeomObject& moveX(double);
   GeomObject& moveY(double);
   GeomObject& moveZ(double);
   // gives subclasses ability to just move

   GeomObject & rotateX(double);
   GeomObject & rotateY(double);
   GeomObject & rotateZ(double);
   // Rotation around the given axis: return pointer to self for convenience

    GeomObject & rotate(double angle, const Vector & axis);
    // Rotation around a specified vector.

    virtual const char*  nameOf()const;
    // makes name available to subclasses using the type info facility

    virtual void printOn(std::ostream& out=std::cout)const;
 
    Command* printCommand(std::ostream* out);
    // create a Command that will have this object print itself, 
    // using the virtual printOn method


protected:
    GeomObject(){}
private:
    class Print;

};

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~inlines~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

inline std::ostream& operator<<(std::ostream& os, const GeomObject& g){
   g.printOn(os); return os;}

#endif

