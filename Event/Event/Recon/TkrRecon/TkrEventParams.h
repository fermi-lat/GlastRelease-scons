#ifndef TkrEventParams_h
#define TkrEventParams_h

#include "geometry/Point.h"
#include "geometry/Vector.h"

#include <iostream>
#include "GaudiKernel/DataObject.h"

#include "GaudiKernel/IInterface.h"

static const CLID& CLID_TkrEventParams = InterfaceID("TkrEventParams", 1, 0);

/** @class TkrEventParams
* @brief Defines the output of the TkrFilterAlg
* 
* 
* $Header$
*/
namespace Event {  // NameSpace

class TkrEventParams : public DataObject 
{    
public:
    
    TkrEventParams() : DataObject(), 
                       m_statusBits(0),
                       m_EventEnergy(0.),
                       m_EventPosition(0.,0.,0.),
                       m_EventAxis(0.,0.,0.) { };

    TkrEventParams(double energy, const Point& pos, const Vector& axis) : DataObject(),
                                                                          m_statusBits(0),
                                                                          m_EventEnergy(energy),
                                                                          m_EventPosition(pos),
                                                                          m_EventAxis(axis) {};

    virtual ~TkrEventParams() { }
    
    /// Retrieve reference to class definition structure
    virtual const CLID& clID() const  { return TkrEventParams::classID(); }
    static const CLID& classID() { return CLID_TkrEventParams; }

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

    //const Point&  getEventPosition() const {return m_EventPosition;}
    //const Vector& getEventAxis()     const {return m_EventAxis;}

    /// Modify data members
    void  setStatusBit(const unsigned int stat) {m_statusBits |= stat;}
    void  setEventEnergy(const double energy)   {m_EventEnergy   = energy;}
    void  setEventPosition(const Point& pos)    {m_EventPosition = pos;}
    void  setEventAxis(const Vector& axis)      {m_EventAxis     = axis;}

private:
    unsigned int m_statusBits;
    double       m_EventEnergy;
    Point        m_EventPosition;
    Vector       m_EventAxis;
};

}

#endif
