/** @file TkrFilterParams.h
*
* $Header$
*
*/

#ifndef TkrFilterParams_h
#define TkrFilterParams_h

#include "geometry/Point.h"
#include "geometry/Vector.h"

#include "GaudiKernel/ObjectList.h"
#include "GaudiKernel/ContainedObject.h"
#include "GaudiKernel/IInterface.h"

#include "Event/RelTable/RelTable.h"

static const CLID& CLID_TkrFilterParams = InterfaceID("TkrFilterParams", 1, 0);

namespace Event {  // NameSpace

/** @class TkrFilterParams
* @brief Defines the output of the TkrFilterAlg
* 
* 
*/

class TkrFilterParams : virtual public ContainedObject 
{    
public:
    
    TkrFilterParams() : m_statusBits(0),
                        m_EventEnergy(0.),
                        m_EventPosition(0.,0.,0.),
                        m_EventAxis(0.,0.,0.),
                        m_numBiLayers(0),
                        m_numIterations(0),
                        m_numHitsTotal(0),
                        m_numDropped(0),
                        m_chiSquare(0.),
                        m_transRms(0.),
                        m_longRmsAve(0.) { };

    TkrFilterParams(double energy, const Point& pos, const Vector& axis) : m_statusBits(0),
                                                                           m_EventEnergy(energy),
                                                                           m_EventPosition(pos),
                                                                           m_EventAxis(axis),
                                                                           m_numBiLayers(0),
                                                                           m_numIterations(0),
                                                                           m_numHitsTotal(0),
                                                                           m_numDropped(0),
                                                                           m_chiSquare(0.),
                                                                           m_transRms(0.),
                                                                           m_longRmsAve(0.) { };

    virtual ~TkrFilterParams() { }
    
    /// Retrieve reference to class definition structure
    virtual const CLID& clID() const  { return TkrFilterParams::classID(); }
    static const CLID& classID() { return CLID_TkrFilterParams; }

    /// Status word bits organized like:
    ///        |  0   0   0   0  |  0   0   0   0  |  0   0   0   0  |  0   0   0   0   |
    ///         <               > <               > <               >  <              >
    enum StatusBits {CALPARAMS  = 0x0001,  //Set if using Calorimeter parameters
                     TKRPARAMS  = 0x0002,  //Set if using Tracker parameters
                     FIRSTPASS  = 0x1000,  //Set if first pass numbers
                     SECONDPASS = 0x2000}; //Set if second pass numbers used

    /// Access data members
    unsigned int getStatusBits()    const {return m_statusBits;}
    double       getEventEnergy()   const {return m_EventEnergy;}
    Point        getEventPosition() const {return m_EventPosition;}
    Vector       getEventAxis()     const {return m_EventAxis;}

    int          getNumBiLayers()   const {return m_numBiLayers;}
    int          getNumIterations() const {return m_numIterations;}
    int          getNumHitsTotal()  const {return m_numHitsTotal;}
    int          getNumDropped()    const {return m_numDropped;}

    double       getChiSquare()     const {return m_chiSquare;}
    double       getTransRms()      const {return m_transRms;}
    double       getLongRmsAve()    const {return m_longRmsAve;}

    /// Modify data members
    void  setStatusBit(const unsigned int stat) {m_statusBits |= stat;}
    void  setEventEnergy(const double energy)   {m_EventEnergy   = energy;}
    void  setEventPosition(const Point& pos)    {m_EventPosition = pos;}
    void  setEventAxis(const Vector& axis)      {m_EventAxis     = axis;}
    void  setNumBiLayers(int numBiLayers)       {m_numBiLayers   = numBiLayers;}
    void  setNumIterations(int numIterations)   {m_numIterations = numIterations;}
    void  setNumHitsTotal(int numHits)          {m_numHitsTotal  = numHits;}
    void  setNumDropped(int numDropped)         {m_numDropped    = numDropped;}
    void  setChiSquare(double chiSquare)        {m_chiSquare     = chiSquare;}
    void  setTransRms(double transRms)          {m_transRms      = transRms;}
    void  setLongRmsAve(double longRmsAve)      {m_longRmsAve    = longRmsAve;}

private:
    unsigned int m_statusBits;
    double       m_EventEnergy;
    Point        m_EventPosition;
    Vector       m_EventAxis;

    // Keep also some diagnostic information
    int          m_numBiLayers;
    int          m_numIterations;
    int          m_numHitsTotal;
    int          m_numDropped;

    // As well as moments output
    double       m_chiSquare;
    double       m_transRms; 
    double       m_longRmsAve;
};

// Typedefs for gaudi container for these objects
typedef ObjectList<TkrFilterParams>                     TkrFilterParamsCol;
typedef TkrFilterParamsCol::iterator                    TkrFilterParamsColPtr;
typedef TkrFilterParamsCol::const_iterator              TkrFilterParamsColConPtr;

class TkrBoundBoxLink;

typedef RelTable<TkrFilterParams, TkrBoundBoxLink>     TkrFilterParamsToLinksTab;
typedef Relation<TkrFilterParams, TkrBoundBoxLink>     TkrFilterParamsToLinksRel;
typedef RelationList<TkrFilterParams, TkrBoundBoxLink> TkrFilterParamsToLinksTabList;

class TkrBoundBox;

typedef RelTable<TkrFilterParams, TkrBoundBox>          TkrFilterParamsToBoxTab;
typedef Relation<TkrFilterParams, TkrBoundBox>          TkrFilterParamsToBoxRel;
typedef RelationList<TkrFilterParams, TkrBoundBox>      TkrFilterParamsToBoxTabList;

class TkrBoundBoxPoint;

typedef RelTable<TkrFilterParams, TkrBoundBoxPoint>     TkrFilterParamsToPointsTab;
typedef Relation<TkrFilterParams, TkrBoundBoxPoint>     TkrFilterParamsToPointsRel;
typedef RelationList<TkrFilterParams, TkrBoundBoxPoint> TkrFilterParamsToPointsTabList;

}

#endif
