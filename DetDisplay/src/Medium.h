/** @file Medium.h
     @brief  declaration of class Medium

   $Header$
*/

#ifndef MEDIUM_H
#define MEDIUM_H

#include "geometry/GeomObject.h"

namespace gui {class DisplayRep; }
class Shape;

/**
*/
class Medium  : public GeomObject
{
public:
    typedef std::string Material;

    /// Standard constructor
    /**
    @param parent The parent: 
    @param vol  The geometric extent of the medium
    @param mat  The material
    @param det  Optional detector, which will be notified by any step inside
    */
    Medium(Medium* parent, Shape* vol, const char* mat= "vacuum");

    /// copy constructor
    Medium(const Medium& old);

    /// destructor
    virtual ~Medium();

    ///  Methods for building Trees of Media. Error if not called from a CompositeMedium
    virtual Medium&  addMedium(Medium* nextMedium);
    ///  Methods for building Trees of Media. Error if not called from a CompositeMedium
    virtual Medium& removeMedium (Medium* oldMedium);


    ///  Methods for setting data members (these all set by constructor and unnessary?)
    Medium&  setParent(Medium* p)       {_parent = p; return *this;}
    Medium&  setVolume(Shape* vol)      {_volume = vol;return *this;}
    Medium&  setTitle(const char *newTitle);

    const std::string& material()const{return _material;}
    ///  Methods for getting data members, attributes
    virtual const char*	nameOf() const;
    virtual int  isComposite()const;
    const char*	title()  const  {return _title? _title : "no title";}

    Shape& 	volume() 	{return *_volume;}
    const Shape& volume()const   {return *_volume;}

    const Medium* getParent() const {return _parent;}
    const Medium* parent()const     {return _parent;}

    ///  co-ordinate transformation
    GeomObject& transform(const CoordTransform& );

    /// pass a gui::DisplayRep object to the Shape for display of the volume
    virtual void createDetectorView(gui::DisplayRep& v){};


protected:
    static unsigned s_count;  // keep informal track of number of Medium's constructed

    // Medium contents:
    Shape*  	_volume;        // Container Volume
    Material 	_material; 	// Material description
    Medium*	_parent;       	// The medium in which this exists, null for root 

    char*  	_title; 	// User setable title

};
//
// Operator overloads for <<
inline std::ostream& operator<< (std::ostream& os, const Medium &m){m.printOn(os);return os;}
inline std::ostream& operator<< (std::ostream& os, const Medium *m){m->printOn(os);return os;}

#endif

