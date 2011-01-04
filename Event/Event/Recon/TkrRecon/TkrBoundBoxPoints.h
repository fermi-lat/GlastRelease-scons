/** @file TkrBoundBoxPoints.h

* @class TkrBoundBoxPoints
*
* @brief This object ties together two TkrVecPoints that are contained within a given bounding box
*        and provides some information about their relationship
*
* last modified 12/08/2010
*
* @authors Tracker Folks
*
* $Header$
*/

#ifndef __TkrBoundBoxPoints_H
#define __TkrBoundBoxPoints_H

#include "Event/Recon/TkrRecon/TkrVecPoint.h"

#include "GaudiKernel/ObjectList.h"
#include "GaudiKernel/ContainedObject.h"
#include "GaudiKernel/IInterface.h"

#include "Event/RelTable/RelTable.h"

// Declare Gaudi object interface ID
static const CLID& CLID_TkrBoundBoxPoint = InterfaceID("TkrBoundBoxPoint",  1, 0);

namespace Event {  // NameSpace

class TkrBoundBoxPoint: virtual public ContainedObject
{
public:
    // constructors
    TkrBoundBoxPoint() : m_parent(0), m_left(0), m_right(0), m_tkrVecPoint(0), m_position(0.,0.,0.) {}
    TkrBoundBoxPoint(const Event::TkrBoundBoxPoint* parent,
                     const Event::TkrBoundBoxPoint* left,
                     const Event::TkrBoundBoxPoint* right,
                     const Event::TkrVecPoint*      tkrVecPoint,
                     const Point&                   position,
                     const double&                  aveSeparation) 
                      : m_parent(parent),
                        m_left(left),
                        m_right(right),
                        m_tkrVecPoint(tkrVecPoint), 
                        m_position(position),
                        m_aveSeparation(aveSeparation)
    {}

    //! Retrieve pointer to class defininition structure
    virtual const CLID& clID()    const   { return TkrBoundBoxPoint::classID(); }
    static  const CLID& classID()         { return CLID_TkrBoundBoxPoint; }

    // destructor
    virtual ~TkrBoundBoxPoint()
    {
        return;
    }
    /// @name Data set methods
    //@{
    /// Set the parent point
    void setBBParent(const Event::TkrBoundBoxPoint* parent)   {m_parent        = parent;}
    /// Set the left daughter
    void setLeft(const Event::TkrBoundBoxPoint* left)         {m_left          = left;}
    /// Set the right daughter
    void setRight(const Event::TkrBoundBoxPoint* right)       {m_right         = right;}
    /// Allow to set the TkrVecPoint
    void setTkrVecPoint(const Event::TkrVecPoint* vecPoint)   {m_tkrVecPoint   = vecPoint;}
    /// Allow to set the position
    void setPosition(const Point& position)                   {m_position      = position;}
    /// Allow to set the average separation up to this point
    void setAveSeparation(const double& aveSeparation)        {m_aveSeparation = aveSeparation;}
    //@}

    /// @name Data access methods
    //@{
    /// Traverse up the tree to the parent
    const Event::TkrBoundBoxPoint* getParent()        const {return m_parent;}
    /// Traverse down to the left daughter
    const Event::TkrBoundBoxPoint* getLeft()          const {return m_left;}
    /// Traverse downn to the right daughter
    const Event::TkrBoundBoxPoint* getRight()         const {return m_right;}
    /// Return the constant pointer to the associated TkrVecPoint (if one)
    const Event::TkrVecPoint*      getTkrVecPoint()   const {return m_tkrVecPoint;}
    /// Return the position of this BB point
    const Point&                   getPosition()      const {return m_position;}
    /// Return the average separation up to this point
    const double                   getAveSeparation() const {return m_aveSeparation;}
    //@}

private:

    // For traversing up and down the binary tree structure
    const Event::TkrBoundBoxPoint* m_parent;
    const Event::TkrBoundBoxPoint* m_left;
    const Event::TkrBoundBoxPoint* m_right;

    // Data members
    const Event::TkrVecPoint*      m_tkrVecPoint;
    Point                          m_position;
    double                         m_aveSeparation;
};    

// Typedefs for gaudi container for these objects
typedef ObjectList<TkrBoundBoxPoint>         TkrBoundBoxPointsCol;
typedef TkrBoundBoxPointsCol::iterator       TkrBoundBoxPointsColPtr;
typedef TkrBoundBoxPointsCol::const_iterator TkrBoundBoxPointsColConPtr;

class TkrBoundBox;

// Typedefs for relating to bounding boxes
typedef RelTable<TkrBoundBoxPoint, TkrBoundBox>     TkrBoundBoxPointsToBoxTab;
typedef Relation<TkrBoundBoxPoint, TkrBoundBox>     TkrBoundBoxPointsToBoxRel;
typedef RelationList<TkrBoundBoxPoint, TkrBoundBox> TkrBoundBoxPointsToBoxTabList;

}; // Namespace

#endif // __TkrBoundBoxPoints_H
