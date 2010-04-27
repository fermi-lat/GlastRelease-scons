#ifndef DiagDataOverlay_h
#define DiagDataOverlay_h

#include <iostream>
#include <vector>

#include "GaudiKernel/Kernel.h"
#include "GaudiKernel/StreamBuffer.h"
#include "GaudiKernel/ContainedObject.h"
#include "GaudiKernel/ObjectVector.h"
#include "GaudiKernel/IInterface.h"

#include "Event/TopLevel/Definitions.h"
#include "Event/TopLevel/EventModel.h"

// From EBF Online library
#include "CALdiagnostic.h"
#include "TKRdiagnostic.h"

/**
* @class DiagDataOverlay
* @brief TDS for storing the diagnostic (trigger primitives)
*        for use in the Overlay
*
* Constains 2 local classes for CAL and TKR trigger primitives
*/

static const CLID& CLID_DiagDataOverlay = InterfaceID("DiagDataOverlay", 1, 0);

namespace Event 
{
    
class CalDiagDataOverlay 
{
public:
    CalDiagDataOverlay(unsigned datum, unsigned tower, unsigned layer) 
        : m_datum(datum), m_tower(tower), m_layer(layer)
    {
    }

    CalDiagDataOverlay()  {}
    ~CalDiagDataOverlay() {}

    unsigned dataWord() const { return m_datum; };
    unsigned tower() const { return m_tower; };
    unsigned layer() const { return m_layer; };

private:
    unsigned m_datum;
    unsigned m_tower;
    unsigned m_layer;
};

class TkrDiagDataOverlay 
{
public:

    TkrDiagDataOverlay(unsigned datum, unsigned tower, unsigned gtcc) 
        : m_datum(datum), m_tower(tower), m_gtcc(gtcc) {}
        TkrDiagDataOverlay()  {}
        ~TkrDiagDataOverlay() {}

        unsigned dataWord() const { return m_datum; };
        unsigned tower() const { return m_tower; }; 
        unsigned gtcc() const { return m_gtcc; };
        // add this access method to make it easier to simulate diagnostic info
        void setDataWord( unsigned int datum) { m_datum = datum; };
private:
    unsigned m_datum;
    unsigned m_tower;  
    unsigned m_gtcc;
};

class DiagDataOverlay : public DataObject
{
public:
    virtual ~DiagDataOverlay();
    virtual std::ostream& fillStream(std::ostream &s) const;
    friend std::ostream& operator << (std::ostream &s, const DiagDataOverlay &obj);

    DiagDataOverlay();

    /// Retrieve reference to class definition structure
    virtual const CLID& clID() const   { return DiagDataOverlay::classID(); }
    static const CLID& classID()       { return CLID_DiagDataOverlay; };

    inline void addCalDiagnostic(const CalDiagDataOverlay& cal) { m_cal.push_back(cal); };
    inline int getNumCalDiagnostic() const { return m_cal.size(); };
    const CalDiagDataOverlay& getCalDiagnosticByIndex(unsigned int ind)const { return m_cal[ind]; };
    inline void addTkrDiagnostic(const TkrDiagDataOverlay& tkr) { m_tkr.push_back(tkr); };
    inline int getNumTkrDiagnostic() const { return m_tkr.size(); };
    const TkrDiagDataOverlay& getTkrDiagnosticByIndex(unsigned int ind)const { return m_tkr[ind]; };

    // needed for simulating diagnostic data
    void setTkrDataWordByIndex( unsigned int ind, unsigned int datum) 
    {
        m_tkr[ind].setDataWord( datum );
    }
private:
    std::vector<CalDiagDataOverlay> m_cal;
    std::vector<TkrDiagDataOverlay> m_tkr;
};

// inline the public get methods for clients besides OnboardFilter.

// the put methods, here for now
inline DiagDataOverlay::DiagDataOverlay()
{
    m_cal.clear();
    m_tkr.clear();
}

inline DiagDataOverlay::~DiagDataOverlay() {}


inline std::ostream& DiagDataOverlay::fillStream(std::ostream &s) const
{
    int i;
    for (i = 0; i < getNumCalDiagnostic(); i++) {
        s<<"Full CAL Data Word # " << i << " " << m_cal[i].dataWord() <<std::endl;
    }
    for (i = 0; i < getNumTkrDiagnostic(); i++) {
        s<<"Full TKR Data Word # " << i << " " <<m_tkr[i].dataWord()<<std::endl;
    }
    return s;
}


inline std::ostream& operator<<(std::ostream &s, const DiagDataOverlay &obj)
{
    return obj.fillStream(s);
}

}//namespace Event 


#endif
