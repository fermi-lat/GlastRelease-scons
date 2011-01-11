/** @file TkrBoundBoxLinks.h

* @class TkrBoundBoxLinks
*
* @brief This object ties together a set of TkrBoundBoxes with their corresponding TkrBoundBoxPoints
*
* last modified 01/10/11
*
* @authors Tracker Folks
*
* $Header$
*/

#ifndef __TkrBoundBoxLinks_H
#define __TkrBoundBoxLinks_H

#include "Event/Recon/TkrRecon/TkrVecPoint.h"

#include "GaudiKernel/ObjectList.h"
#include "GaudiKernel/ContainedObject.h"
#include "GaudiKernel/IInterface.h"

// Declare Gaudi object interface ID
static const CLID& CLID_TkrBoundBoxLink = InterfaceID("TkrBoundBoxLink",  1, 0);

namespace Event {  // NameSpace

class TkrBoundBoxPoint;
class TkrBoundBox;

class TkrBoundBoxLink : virtual public ContainedObject
{
public:
    // constructors
    TkrBoundBoxLink() : m_parentLink(0), m_topPoint(0), m_associatedBox(0) {}
    TkrBoundBoxLink(const Event::TkrBoundBoxLink*  parent,
                    const Event::TkrBoundBoxPoint* topPoint,
                    const Event::TkrBoundBox*      boundBox,
                    const Point&                   position,
                    const double&                  distToParent)
                     : m_parentLink(parent),
                       m_topPoint(topPoint),
                       m_associatedBox(boundBox),
                       m_position(position),
                       m_distToParent(distToParent)
    {}

    //! Retrieve pointer to class defininition structure
    virtual const CLID& clID()    const   { return TkrBoundBoxLink::classID(); }
    static  const CLID& classID()         { return CLID_TkrBoundBoxLink; }

    // destructor
    virtual ~TkrBoundBoxLink()
    {
        return;
    }
    /// @name Data set methods
    //@{
    /// Set the parent point
    void setParentLink(const Event::TkrBoundBoxLink* parent)  {m_parentLink    = parent;}
    /// Set the top point associated to this link
    void setTopPoint(const Event::TkrBoundBoxPoint* topPoint) {m_topPoint      = topPoint;}
    /// Set the associated TkrBoundBox
    void setBoundBox(const Event::TkrBoundBox* boundBox)      {m_associatedBox = boundBox;}
    /// Allow to set the position
    void setPosition(const Point& position)                   {m_position      = position;}
    /// Allow to set the average separation up to this point
    void setDistToParent(const double& distToParent)          {m_distToParent  = distToParent;}
    //@}

    /// @name Data access methods
    //@{
    /// Return the parent for this link
    const Event::TkrBoundBoxLink*  getParent()        const {return m_parentLink;}
    /// Return the top point associated to this link
    const Event::TkrBoundBoxPoint* getTopPoint()      const {return m_topPoint;}
    /// Return the TkrBoundBox associated to this link
    const Event::TkrBoundBox*      getBoundBox()      const {return m_associatedBox;}
    /// Return the position of this link (the box position)
    const Point&                   getPosition()      const {return m_position;}
    /// Return the distance to the parent link
    const double                   getDistToParent()  const {return m_distToParent;}
    //@}

private:

    // For traversing up and down the binary tree structure
    const Event::TkrBoundBoxLink*  m_parentLink;
    const Event::TkrBoundBoxPoint* m_topPoint;
    const Event::TkrBoundBox*      m_associatedBox;
    Point                          m_position;
    double                         m_distToParent;
};    

// Typedefs for gaudi container for these objects
typedef ObjectList<TkrBoundBoxLink>         TkrBoundBoxLinksCol;
typedef TkrBoundBoxLinksCol::iterator       TkrBoundBoxLinksColPtr;
typedef TkrBoundBoxLinksCol::const_iterator TkrBoundBoxLinksColConPtr;

class TkrBoundBox;

#include "Event/RelTable/RelTable.h"

// Typedefs for relating to bounding boxes
typedef RelTable<TkrBoundBoxLink,     TkrBoundBox> TkrBoundBoxLinksToBoxTab;
typedef Relation<TkrBoundBoxLink,     TkrBoundBox> TkrBoundBoxLinksToBoxRel;
typedef RelationList<TkrBoundBoxLink, TkrBoundBox> TkrBoundBoxLinksToBoxTabList;

}; // Namespace

#endif // __TkrBoundBoxLinks_H
