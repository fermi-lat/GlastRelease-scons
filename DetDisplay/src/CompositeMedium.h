//  $Id$
//
// This file is part of Gismo 2
//
#ifndef COMPOSITEMEDIUM_H
#define COMPOSITEMEDIUM_H


#include "Medium.h"

#include <vector>

typedef std::vector<Medium*> MediumList;

/**
 Inherits from Medium to:
<UL>
   <LI> define behavior for Mediums having internal Mediums, often
	simply passing the same message
   <LI> store the internal Mediums
   <LI> implement child-related operations addMedium, removeMedium and
	provide for iteration over the list
</UL>
*/
class CompositeMedium  : public Medium , public MediumList
{
 public:

    ///  Default Constructor 
   CompositeMedium(Medium* =0, float size=100);

   /// Standard constructor
   /**
      @param mother The mother medium, which must be composite. This medium will be added to its children. The volume
      @param vol    Volume assigned to this medium
      @param mat    Material
      @param det    OPtional detector object, that will be notified for each step in this Medium
    */
   CompositeMedium(Medium* mother, Shape* vol, const char* mat= "vacuum", Detector * det= 0);
//// destructor
   virtual ~CompositeMedium();

///   child-related operations
   virtual Medium&  addMedium(Medium* nextMedium);
///   child-related operations
   virtual Medium& removeMedium (Medium* oldMedium);
///   child-related operations
   void deleteInnerMedia();

///   child-related operations
   int innerMediaCount()const  {return size();}
///   child-related operations
   unsigned int size()const { return MediumList::size();}

///   child-related operations
   const Medium* innerMedium(int i)const{return operator[](i);}
///   child-related operations
   const char* nameOf() const;
///   child-related operations
   int isComposite() const;
#if 0
   ///  set tracking attributes:
   virtual Medium&	setKECutOff(float keCut);
   ///  set tracking attributes:
   virtual Medium&	setMaxStep(float mxStep);
   ///  set tracking attributes:
   virtual Medium&	setField(Field* );
#endif
    ///   co-ordinate transformation

   GeomObject& transform(const CoordTransform& );
#if 0
///   Methods for tracking
   virtual double distanceToLeave( const Ray& r, 
				   const Medium*& 
				   newstuff,
				   double maxStep ) const;
   /// will check if the Ray intersects an internal Medium

   /// return the Medium that the point is inside, perhaps an internal one
   virtual const Medium * inside(const Point& r)const;


   virtual void clear();
   virtual void generateResponse();
   virtual void accept(DetectorVisitor&);
   virtual void readData(std::istream&);
   virtual void writeData(std::ostream&);
   virtual void notify();


///     Methods associated with GUI and Printing and I/O
   virtual void printOn( std::ostream& os = std::cout ) const;
///     Methods associated with GUI and Printing and I/O
   virtual void printResponse(std::ostream& = std::cout) const;
#endif
//     Methods associated with GUI and Printing and I/O

   /// pass a gui::DisplayRep object to the Shape  and invoke for all inner media
   virtual void createDetectorView(gui::DisplayRep& v);
#if 0 
   /// pass a gui::DisplayRep object to the Detector and invoke for all inner media
   virtual void createResponseView(gui::DisplayRep& v);
#endif

 private:

};


#endif

