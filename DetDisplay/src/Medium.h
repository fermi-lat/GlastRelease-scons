// $Header$
//


#ifndef GISMO_MEDIUM_H
#define GISMO_MEDIUM_H

#ifdef __GNUG__
#pragma interface
#endif

#include "geometry/GeomObject.h"

class Material;
class Detector;
namespace gui {class DisplayRep; }
class Field;
class Shape;
class Ray;
class Vector;
class Point;
class DetectorVisitor;

/**
 Medium is an interface class for MCParticle, combining information on the material, field, and geometry.  
<BR>
  This information is use in the
propagate method of  MCParticle. Furthermore mediums provide "hooks" to 
generate detector responses and do event anaylsis.
<BR>
  Mediums contain mediums and hence form a hierarchical arrangement of the 
  geometry.  In particular a  CompositeMedium has a list of mediums it contains
  as well as being in the list of mediums of its parent. 
*/
class Medium  : public GeomObject
{
public:
    
    /// default  Constructors / destructor
    Medium(Medium* parent =0, float size=100);
    /// Standard constructor
    /**
    @param parent The parent: insure integrity of the medium hierarchy.
    and to make child tracking attributes default to those of the parent
    @param vol  The geometric extent of the medium
    @param mat  The material
    @param det  Optional detector, which will be notified by any step inside
    */
    Medium(Medium* parent, Shape* vol, const char* mat= "vacuum", Detector * det= 0);

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

    
    
    ///  Methods for getting data members, attributes
    virtual const char*	nameOf() const;
    virtual int  isComposite()const;
    const char*	title()  const  {return _title? _title : "no title";}
    
    Shape& 	volume() 	{return *_volume;}
    const Shape& volume()const   {return *_volume;}
    Field& 	field()const    {return *_field;}
    Material& 	material()const {return *_material;}
    
    const Medium* getParent() const {return _parent;}
    const Medium* parent()const     {return _parent;}
    
    float 	kECutOff()const  {return _keCutOff;}
    float 	maxStep()const   {return _maxStep;}
    Detector* 	detector()const  {return _detector;}
    
    ///  co-ordinate transformation
        GeomObject& transform(const CoordTransform& );
    

    /// pass a gui::DisplayRep object to the Shape for display of the volume
    virtual void createDetectorView(gui::DisplayRep& v);

    
protected:
    void set_defaults();
    static unsigned s_count;  // keep informal track of number of Medium's constructed

    // Medium contents:

    Shape*  	_volume;        // Container Volume
    Material* 	_material; 	// Material description
    Detector*	_detector;	// (optional) Detector
    Field*  	_field;         // Magnetic Field, null for no field
    
    Medium*	_parent;       	// The medium in which this exists, null for root 
    
    // attributes for tracking
    float 	_keCutOff;	// Cut-off momentum for swimming
    float 	_maxStep;	// Step size limit
    
    char*  	_title; 	// User setable title
    
    static const Medium *lastMedium; // pointer to the last distanceToLeave Medium
       
};
//
// Operator overloads for <<
inline std::ostream& operator<< (std::ostream& os, const Medium &m){m.printOn(os);return os;}
inline std::ostream& operator<< (std::ostream& os, const Medium *m){m->printOn(os);return os;}

#endif

