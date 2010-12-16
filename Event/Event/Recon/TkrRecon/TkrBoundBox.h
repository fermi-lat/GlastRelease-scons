/** @file TkrBoundBox.h

* @class TkrBoundBox
*
* @brief This object defines a bounding box around a group of TkrVecPoints in a given bilayer. 
*        The point of this is to facilitate the determination of the main event axis in the tracker
*
* last modified 12/08/2010
*
* @authors Tracker Folks
*
* $Header$
*/

#ifndef __TkrBoundBox_H
#define __TkrBoundBox_H

#include "Event/Recon/TkrRecon/TkrVecPoint.h"
#include <vector>

#include "GaudiKernel/ObjectList.h"
#include "GaudiKernel/ContainedObject.h"
#include "GaudiKernel/IInterface.h"

#include "geometry/Point.h"

// Declare Gaudi object interface ID
static const CLID& CLID_TkrBoundBox = InterfaceID("TkrBoundBox",  1, 0);

namespace Event {  // NameSpace

typedef std::list<const TkrVecPoint*> TkrVecPointConstList;

class TkrBoundBox: public TkrVecPointConstList, virtual public ContainedObject
{
public:
    // constructors
    TkrBoundBox()
    {}

    //! Retrieve pointer to class defininition structure
    virtual const CLID& clID()    const   { return TkrBoundBox::classID(); }
    static  const CLID& classID()         { return CLID_TkrBoundBox; }

    // destructor
    virtual ~TkrBoundBox()
    {
        return;
    }
    /// @name Data set methods
    //@{
    /// Set the biLayer for this box
    void setBiLayer(const int& biLayer)       {m_biLayer = biLayer;}
    /// Set the low corner for this box
    void setLowCorner(const Point& pt)        {m_lowCorner = pt;}
    /// Set the high corner for this box
    void setHighCorner(const Point& pt)       {m_highCorner = pt;}
    /// Set the average position of the TkrVecPoints in the box
    void setAveragePosition(const Point& pt)  {m_averagePos = pt;}
    /// Set the hit density
    void setHitDensity(const double& density) {m_hitDensity = density;}
    /// Set the mean distance between TkrVecPoints
    void setMeanDist(const double& meanDist)  {m_meanDist = meanDist;}
    /// Set the rms distance between TkrVecPoints
    void setRmsDist(const double& rmsDist)    {m_rmsDist = rmsDist;}
    //@}

    /// @name Data access methods
    //@{
    /// Return the biLayer this box belongs too
    const int     getBiLayer()         const {return m_biLayer;}
    /// Get the low corner of the bounding box
    const Point&  getLowCorner()       const {return m_lowCorner;}
    /// Get the high corner of the bounding box
    const Point&  getHighCorner()      const {return m_highCorner;}
    /// Average position of the TkrVecPoints contained in the box
    const Point&  getAveragePosition() const {return m_averagePos;}
    /// Center of the bounding box
    const Point&  getBoxCenterPos();
    /// Return the "hit density" of TkrVecPoints for this box
    const double  getHitDensity()      const {return m_hitDensity;}
    /// Return the mean distance between TkrVecPoints from MST
    const double  getMeanDist()        const {return m_meanDist;}
    /// Return the rms distance between TkrVecPoints from MST
    const double  getRmsDist()         const {return m_rmsDist;}
    /// Return number of TkrVecPoints in the box
    const int     getNumPoints()       const {return size();}
    //@}

private:

    // data members
    int    m_biLayer;          // The biLayer this box belongs too
    
    Point  m_lowCorner;        // The left bottom corner (in a positive system)
    Point  m_highCorner;       // The right upper corner

    Point  m_averagePos;       // The average position of the contained TkrVecPoints
    Point  m_boxCenterPos;     // The coordinates of the center of the box

    double m_hitDensity;       // The "density" of the contained hits

    double m_meanDist;         // The mean distance from the MST for TkrVecPoints in this box
    double m_rmsDist;          // The rms distance from the MST for TkrVecPoints in this box
};

const Point& TkrBoundBox::getBoxCenterPos()
{
    m_boxCenterPos  = m_lowCorner;
    m_boxCenterPos += m_highCorner;
    m_boxCenterPos *= 0.5;

    return m_boxCenterPos;
}
    

// Typedefs for gaudi container for these objects
typedef ObjectList<TkrBoundBox>        TkrBoundBoxCol;
typedef TkrBoundBoxCol::iterator       TkrBoundBoxColPtr;
typedef TkrBoundBoxCol::const_iterator TkrBoundBoxColConPtr;

}; // Namespace

#endif // __TkrBoundBox_H
