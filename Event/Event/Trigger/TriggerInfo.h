#ifndef EVENT_TriggerInfo_H
#define EVENT_TriggerInfo_H

#include "GaudiKernel/Kernel.h"
#include "GaudiKernel/StreamBuffer.h"
#include "GaudiKernel/ContainedObject.h"
#include "GaudiKernel/ObjectVector.h"

#include "Event/TopLevel/Definitions.h"
#include "Event/TopLevel/EventModel.h"

#include <iostream>
#include <map>

/**
* @class TriggerInfo
* @brief TDS for TriggerInfo data
*/
static const CLID& CLID_TriggerInfo = InterfaceID("TriggerInfo", 1, 0);

namespace Event {

/** @class TriggerInfo
  * @brief Local storage of TriggerInfo data
  * $Header$
*/
class TriggerInfo : public DataObject
{
public:
    /// Define the TileList map
    typedef std::map<std::string, unsigned int> TileList;

    TriggerInfo() { clear(); };
    virtual ~TriggerInfo() { clear(); };

    /// Retrieve reference to class definition structure
    virtual const CLID& clID()    const { return TriggerInfo::classID(); }
    static  const CLID& classID()       { return CLID_TriggerInfo; }

    void clear();

    /// Global initialization method
    void initTriggerInfo(unsigned int    m_triggerBits,
                         unsigned short  tkrVector,
                         unsigned short  roi,
                         unsigned short  calLE,
                         unsigned short  calHE,
                         unsigned short  cno,
                         unsigned short  deltaEventTime,
                         unsigned short  deltaWindowTime,
                         const TileList& tileList);

    /// Setindividual members
    void setTriggerBits(unsigned int triggerBits)    {m_triggerBits         = triggerBits;}
    void setTkrVector(unsigned short tkrVector)      {m_tkrVector           = tkrVector;}
    void setRoiVector(unsigned short roi)            {m_roiVector           = roi;}
    void setCalLeVector(unsigned short calLE)        {m_calLeVector         = calLE;}
    void setCalHeVector(unsigned short calHE)        {m_calHeVector         = calHE;}
    void setCnoVector(unsigned short cno)            {m_cnoVector           = cno;}
    void setDeltaEventTime(unsigned short time)      {m_deltaEventTime      = time;}
    void setDeltaWindowOpenTime(unsigned short time) {m_deltaWindowOpenTime = time;}
    void setTileList(const TileList& tileList)       {m_tileList            = tileList;}

    /// Retrieve values
    unsigned int    getTriggerBits()         const {return m_triggerBits;}
    unsigned short  getTkrVector()           const {return m_tkrVector;}
    unsigned short  getRoiVector()           const {return m_roiVector;}
    unsigned short  getCalLeVector()         const {return m_calLeVector;}
    unsigned short  getCalHeVector()         const {return m_calHeVector;}
    unsigned short  getCnoVector()           const {return m_cnoVector;}
    unsigned short  getDeltaEventTime()      const {return m_deltaEventTime;}
    unsigned short  getDeltaWindowOpenTime() const {return m_deltaWindowOpenTime;}
    const TileList& getTileList()            const {return m_tileList;}

    virtual std::ostream& fillStream(std::ostream &s) const;
    friend std::ostream& operator << (std::ostream &s, const TriggerInfo &obj);

private:

    unsigned int   m_triggerBits;
    unsigned short m_tkrVector;
    unsigned short m_roiVector;
    unsigned short m_calLeVector;
    unsigned short m_calHeVector;
    unsigned short m_cnoVector;
    unsigned short m_deltaEventTime;
    unsigned short m_deltaWindowOpenTime;
    TileList       m_tileList;
};


inline void TriggerInfo::initTriggerInfo(unsigned int    triggerBits,
                                         unsigned short  tkrVector,
                                         unsigned short  roi,
                                         unsigned short  calLE,
                                         unsigned short  calHE,
                                         unsigned short  cno,
                                         unsigned short  deltaEventTime,
                                         unsigned short  deltaWindowTime,
                                         const TileList& tileList)
{
    m_triggerBits         = triggerBits;
    m_tkrVector           = tkrVector;
    m_roiVector           = roi;
    m_calLeVector         = calLE;
    m_calHeVector         = calHE;
    m_cnoVector           = cno;
    m_deltaEventTime      = deltaEventTime;
    m_deltaWindowOpenTime = deltaWindowTime;

    m_tileList            = tileList;

    return;
}

inline void TriggerInfo::clear() 
{
    m_triggerBits         = 0;
    m_tkrVector           = 0;
    m_roiVector           = 0;
    m_calLeVector         = 0;
    m_calHeVector         = 0;
    m_cnoVector           = 0;
    m_deltaEventTime      = 0xffff;
    m_deltaWindowOpenTime = 0xffff;

    m_tileList.clear();

    return;
}



inline std::ostream& TriggerInfo::fillStream(std::ostream &s) const
{
    s << "TriggerInfo:" <<std::endl;
/*
    s << "ROI vector = 0x" << std::hex << std::setw(4) << std::setfill('0') << roiVector() << std::endl;
    s << "TKR vector = 0x" << std::hex << std::setw(4) << tkrVector() << std::endl;
    s << "CAL HE vector = 0x" << std::hex << std::setw(4) << m_cal_HE_Vector << std::endl;
    s << "CAL LE vector = 0x" << std::hex << std::setw(4) << m_cal_LE_Vector << std::endl;
    s << "Condition Summary = 0x" << std::hex << std::setw(4) << m_conditionSummary << std::endl;
    s << "CNO vector        = 0x" << std::hex << std::setw(4) << m_cno_Vector << std::endl;
    m_tileList.fillStream(s);
    s << "Live time         = 0x" << std::hex << std::setw(8) << m_liveTime << std::dec << " = "; 
    s << m_liveTime << std::endl;
    s << "Prescaled         = 0x" << std::hex << std::setw(8) << m_prescaled << std::dec << " = "; 
    s << m_prescaled << std::endl;
    s << "Discarded         = 0x" << std::hex << std::setw(8) << m_discarded << std::dec << " = ";
    s << m_discarded << std::endl;
    //s << "Sent              = 0x" << std::hex << std::setw(8) << m_sent << std::dec << " = " ;
    //s << m_sent << std::endl;
    s << "Trigger Time      = 0x" << std::hex << std::setw(8) << m_triggerTime << std::dec << " = " ;
    s << m_triggerTime << std::endl;
    m_onePpsTime.fillStream(s);
    s << "Delta event time = 0x" << std::hex << std::setw(8) << m_deltaEventTime << std::dec << " = ";
    s << m_deltaEventTime << std::endl;
*/
    return s;
}

inline std::ostream& operator<<(std::ostream &s, const TriggerInfo &obj)
{
    return obj.fillStream(s);
}

}//namespace Event


#endif
