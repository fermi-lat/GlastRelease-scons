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
static const CLID& CLID_TkrBoundBoxPoints = InterfaceID("TkrBoundBoxPoints",  1, 0);

namespace Event {  // NameSpace

class TkrBoundBoxPoints: virtual public ContainedObject
{
public:
    // constructors
    TkrBoundBoxPoints() : m_firstVecPoint(0), m_secondVecPoint(0), m_separation(0.) {}
    TkrBoundBoxPoints(const Event::TkrVecPoint* firstVecPoint,
                      const Event::TkrVecPoint* secondVecPoint,
                      const double              separation) 
                      : m_firstVecPoint(firstVecPoint), 
                        m_secondVecPoint(secondVecPoint), 
                        m_separation(separation) 
    {}

    //! Retrieve pointer to class defininition structure
    virtual const CLID& clID()    const   { return TkrBoundBoxPoints::classID(); }
    static  const CLID& classID()         { return CLID_TkrBoundBoxPoints; }

    // destructor
    virtual ~TkrBoundBoxPoints()
    {
        return;
    }
    /// @name Data set methods
    //@{
    /// Allow to set the first TkrVecPoint
    void setFirstTkrVecPoint(const Event::TkrVecPoint* vecPoint)  {m_firstVecPoint = vecPoint;}
    /// Allow to set the second TkrVecPoint
    void setSecondTkrVecPoint(const Event::TkrVecPoint* vecPoint) {m_firstVecPoint = vecPoint;}
    /// Allow to set the separation between them
    void setSeparation(const double separation)                   {m_separation    = separation;}
    //@}

    /// @name Data access methods
    //@{
    /// Return the constant pointer to the first TkrVecPoint
    const Event::TkrVecPoint* getFirstTkrVecPoint()  const {return m_firstVecPoint;}
    /// Return the constant pointer to the second TkrVecPoint
    const Event::TkrVecPoint* getSecondTkrVecPoint() const {return m_secondVecPoint;}
    /// Return the separation between the two
    const double              getSeparation()        const {return m_separation;}
    //@}

private:

    // data members
    const Event::TkrVecPoint* m_firstVecPoint;
    const Event::TkrVecPoint* m_secondVecPoint;
    double                    m_separation;
};    

// Typedefs for gaudi container for these objects
typedef ObjectList<TkrBoundBoxPoints>        TkrBoundBoxPointsCol;
typedef TkrBoundBoxPointsCol::iterator       TkrBoundBoxPointsColPtr;
typedef TkrBoundBoxPointsCol::const_iterator TkrBoundBoxPointsColConPtr;

class TkrBoundBox;

typedef RelTable<TkrBoundBoxPoints, TkrBoundBox>     TkrBoundBoxPointsToBoxTab;
typedef Relation<TkrBoundBoxPoints, TkrBoundBox>     TkrBoundBoxPointsToBoxRel;
typedef RelationList<TkrBoundBoxPoints, TkrBoundBox> TkrBoundBoxPointsToBoxTabList;

}; // Namespace

#endif // __TkrBoundBoxPoints_H
