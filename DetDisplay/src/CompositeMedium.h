//  $Id$
//
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

    /// Standard constructor
    /**
    @param mother The mother medium, which must be composite. This medium will be added to its children. The volume
    @param vol    Volume assigned to this medium
    @param mat    Material
    @param det    OPtional detector object, that will be notified for each step in this Medium
    */
    CompositeMedium(Medium* mother, Shape* vol, const char* mat= "vacuum");
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

    ///   co-ordinate transformation

    GeomObject& transform(const CoordTransform& );
    //     Methods associated with GUI and Printing and I/O

    /// pass a gui::DisplayRep object to the Shape  and invoke for all inner media
    virtual void createDetectorView(gui::DisplayRep& v);

};


#endif

