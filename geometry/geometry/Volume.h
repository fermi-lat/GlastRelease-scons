// $Id$
//
//

#ifndef __VOLUME_H
#define __VOLUME_H

#include "geometry/Shape.h"
#include "geometry/Surface.h"
#include <string> //gcc bug: shield vector by string
#include <vector>
typedef std::vector<Surface*> SurfaceList;

/** Defines subclass of Shape that uses a list of Surfaces to define the
 the volume. This class implements methods associated with the surfaces
 */
class Volume : public  Shape , protected SurfaceList
{
 protected:
   explicit Volume(unsigned n);
   // default constructor (with number of surfaces to allocate), destructor

 public:
      ~Volume();

   /** return the distance along the Ray to leave (if inside ), or
       enter (if outside). Can be negative if the origin of the Ray is
       close to the surface. testDist is the maxiumum distance required.
       return MAX_FLT if no intersection between the Ray and the volume
   */
   virtual double distanceToLeave( const Ray& r, double testDist)const;
   virtual double distanceToEnter( const Ray& r, double testDist)const;

   /// is the point inside the volume?
   virtual int inside ( const Point& x ) const;

   /// rotate, translate
   GeomObject& transform(const CoordTransform&);

   virtual void printOn( std::ostream& os = std::cout ) const;

   const Surface& surface(int i)const { return *operator[](i); }
   unsigned surfaceCount()const{return size();}
   // public access to surfaces


    virtual int getBoundaryIndex(const Point& p)const;
    // allows client using distanceToLeave to associate an index with a boundary,
    // identified by the point p

    virtual int lastBoundaryIndex()const;
    // distanceToLeave sets index of left boundary: client can use

    double getMaxDimension() const { return max_dimension; }


 protected:
   // following methods accessible to subclasses only
   Surface& surface(int i) { return *operator[](i); }
   void addSurface(Surface* s);
   void deleteSurfaces();

   virtual void calcMaxDimension() { max_dimension = FLT_MAX; }

   // protected data members
   double max_dimension;	// maximum dimension of a Volume

 private:

};

//
// ---------------------- Inlines ----------------------------------
inline std::ostream& operator<<( std::ostream& os, const Volume& v )
 {v.printOn(os); return os; }
inline std::ostream& operator<<( std::ostream& os, const Volume* v )
 {v->printOn(os); return os;}

#endif

