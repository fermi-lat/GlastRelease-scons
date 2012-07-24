#ifndef GemOverlay_h
#define GemOverlay_h

#include <iostream>
#include <vector>

#include "GaudiKernel/Kernel.h"
#include "GaudiKernel/StreamBuffer.h"
#include "GaudiKernel/ContainedObject.h"
#include "GaudiKernel/ObjectVector.h"
#include "GaudiKernel/IInterface.h"

/**
* @class GemOverlay
* @brief Overlay TDS for GEM data
*
* $Header$
*/
static const CLID& CLID_GemOverlay = InterfaceID("GemOverlay", 1, 0);

namespace Event {

class GemOverlayTileList
{
public:
  GemOverlayTileList()  {}
  GemOverlayTileList (unsigned short xzm, unsigned short xzp, unsigned short yzm,
                      unsigned short yzp, unsigned xy, unsigned short rbn, 
                      unsigned short na) :
                      m_xzp(xzp), m_xzm(xzm), m_yzp(yzp), m_yzm(yzm), m_xy(xy),
                      m_rbn(rbn), m_na(na) { };

  GemOverlayTileList(const GemOverlayTileList &tileList) : 
                     m_xzp(tileList.m_xzp), m_xzm(tileList.m_xzm), m_yzp(tileList.m_yzp), 
                     m_yzm(tileList.m_yzm), m_xy(tileList.m_xy), m_rbn(tileList.m_rbn), 
                     m_na(tileList.m_na) { };
               
  ~GemOverlayTileList() {}

  void init (unsigned short xzm, unsigned short xzp, unsigned short yzm,
               unsigned short yzp, unsigned xy, unsigned short rbn, 
               unsigned short na) { 
               m_xzp = xzp;
               m_xzm = xzm;
               m_yzp = yzp;
               m_yzm = yzm;
               m_xy = xy;
               m_rbn = rbn;
               m_na = na;
               };

  void clear() {
               m_xzp = 0;
               m_xzm = 0;
               m_yzp = 0;
               m_yzm = 0;
               m_xy = 0;
               m_rbn = 0;
               m_na = 0;
  };

    inline std::ostream& fillStream(std::ostream &s) const{
        s << "tile list: " << std::endl;
		s << "XZM, XZP = 0x" << std::hex << std::setw(4) << getXzm() << ", 0x"; 
        s << std::hex << getXzp() << std::endl;
        s << "YZM, YZP = " << getYzm() << ", " << getYzp() << std::endl;
		s << "XY = " << std::setw(8) << getXy() << std::endl;
		s << "RBN, NA = " << std::setw(4) << getRbn() << ", " << getNa() << std::endl;
        return s;
    }

  unsigned short getXzm()  const { return m_xzm; }
  unsigned short getXzp()  const { return m_xzp; }
  unsigned short getYzm()  const { return m_yzm; };
  unsigned short getYzp()  const { return m_yzp; };
  unsigned       getXy()   const { return m_xy; };
  unsigned short getRbn()  const { return m_rbn; };
  unsigned short getNa()   const { return m_na; };

private:
  unsigned short m_xzp;
  unsigned short m_xzm;
  unsigned short m_yzp;
  unsigned short m_yzm;
  unsigned       m_xy;
  unsigned short m_rbn;
  unsigned short m_na;
};


/** @class GemOverlayOnePpsTime
 *  @brief This is a helper class for parsing the One PPS time structure.
 */
class GemOverlayOnePpsTime
{
public:
  GemOverlayOnePpsTime()  {}
  GemOverlayOnePpsTime(unsigned time, unsigned sec) { 
     m_timebase = time; 
     m_seconds = sec;
  };
  GemOverlayOnePpsTime(const GemOverlayOnePpsTime &time) : 
              m_timebase(time.m_timebase), m_seconds(time.m_seconds) {};

  ~GemOverlayOnePpsTime() {}

  void init(unsigned timebase, unsigned sec) { 
       m_timebase = timebase; 
       m_seconds = sec; 
  };

  void clear() { 
     m_timebase = 0; 
     m_seconds = 0; 
  };

  inline std::ostream& fillStream(std::ostream &s) const{
	  s << "OnePPSTime:" << std::endl;
	  s << "seconds, timebase = " << m_seconds << ", " << m_timebase << std::endl;
      return s;
  }

  unsigned getTimebase()  const { return m_timebase; };
  unsigned getSeconds()   const { return m_seconds; };

private:
  unsigned m_timebase;
  unsigned m_seconds;
};


class GemOverlayDataCondArrivalTime 
{
public:

  void init(unsigned val) { m_condArr = val; };
  unsigned short external()  const;
  unsigned short cno()       const;
  unsigned short calLE()     const;
  unsigned short calHE()     const;
  unsigned short tkr()       const;
  unsigned short roi()       const;
  unsigned condArr() const { return m_condArr; };

private:
  unsigned m_condArr;
  
};




/** @class Gem
  * @brief Local storage of GEM data
  * $Header$
*/
class GemOverlay : public DataObject
{
public:

    /// Defines the 8 bits in the Condition Summary word
    /// Please see Section 1.7.2 The Condition Summary in the 
    /// "The GLT Electronics Module" available at:
    /// http://www-glast.slac.stanford.edu/IntegrationTest/ONLINE/docs/GEM.pdf
    typedef enum {
        ROI = 1,
        TKR = 2,
        CALLE = 4,
        CALHE = 8,
        CNO = 16,
        PERIODIC = 32,
        SOLICITED = 64,
        EXTERNAL = 128 
    } Summary;

    GemOverlay() { clear(); };
    virtual ~GemOverlay() { clear(); };
    virtual std::ostream& fillStream(std::ostream &s) const;
    friend std::ostream& operator << (std::ostream &s, const GemOverlay &obj);


   // GEM doc v2.5 redefined sent() field as condition_arrival_times
    // We keep both accessors to this field for backward compatibility,
    // but overload the memory location with a union.
    union sent_condArr {
        unsigned sent;
        GemOverlayDataCondArrivalTime condArr;
    };

    /// Retrieve reference to class definition structure
    virtual const CLID& clID() const   { return GemOverlay::classID(); }
    static const CLID& classID()       { return CLID_GemOverlay; }

    void clear();

    void initTrigger(unsigned short tkr, unsigned short roi, 
                     unsigned short calLE,
                     unsigned short calHE, unsigned short cno, 
                     unsigned short conditionSummary, 
                     unsigned short missed,
                     const GemOverlayTileList &tileList);

    void setTkrVector( unsigned short tkr);

    void initSummary(unsigned liveTime, unsigned prescaled, 
                     unsigned discarded, GemOverlayDataCondArrivalTime condArr,
                     unsigned triggerTime, const GemOverlayOnePpsTime &time, 
                     unsigned short deltaEvtTime, unsigned short deltaWindOpenTime);
/*
        void initSummary(unsigned liveTime, unsigned prescaled, 
                         unsigned discarded, unsigned sent, 
                         unsigned triggerTime, const GemOverlayOnePpsTime &time, 
                         unsigned short deltaEvtTime, unsigned short deltaWindOpenTime); */


     unsigned short                getTkrVector()           const { return m_tkrVector;};
     unsigned short                getRoiVector()           const { return m_roiVector;};
     unsigned short                getCalLEvector()         const { return m_cal_LE_Vector;};
     unsigned short                getCalHEvector()         const { return m_cal_HE_Vector;};
     unsigned short                getCnoVector()           const { return m_cno_Vector;};
     unsigned short                getConditionSummary()    const { return m_conditionSummary;};
     unsigned short                getMissed()              const { return m_missed; };
     const GemOverlayTileList&     getTileList()            const { return m_tileList; };
     unsigned                      getLiveTime()            const { return m_liveTime; };
     unsigned                      getPrescaled()           const { return m_prescaled;};
     unsigned                      getDiscarded()           const { return m_discarded;};
     unsigned                      getSent()                const { return m_sent_condArr.sent; };
     GemOverlayDataCondArrivalTime getCondArrTime()         const { return m_sent_condArr.condArr; };
     unsigned                      getTriggerTime()         const { return m_triggerTime; };
     const GemOverlayOnePpsTime&   getOnePpsTime()          const { return m_onePpsTime; };
     unsigned short                getDeltaEventTime()      const { return m_deltaEventTime; };
     unsigned short                getDeltaWindowOpenTime() const { return m_deltaWindOpenTime; };

     /// Methods to query bits in the condition summary word
     bool roiSet() const { return( (m_conditionSummary & ROI) != 0); };
     bool tkrSet() const { return( (m_conditionSummary & TKR) != 0); };
     bool calLeSet() const { return( (m_conditionSummary & CALLE) != 0); };
     bool calHeSet() const { return( (m_conditionSummary & CALHE) != 0); };
     bool cnoSet() const { return( (m_conditionSummary & CNO) != 0); };
     bool periodicSet() const { return( (m_conditionSummary & PERIODIC) != 0); };
     bool solicitedSet() const { return( (m_conditionSummary & SOLICITED) != 0) ; };
     bool externalSet() const { return( (m_conditionSummary & EXTERNAL) != 0); };

 private:

     unsigned short       m_tkrVector;
     unsigned short       m_roiVector;
     unsigned short       m_cal_LE_Vector;
     unsigned short       m_cal_HE_Vector;
     unsigned short       m_cno_Vector;
     unsigned short       m_conditionSummary;
     GemOverlayTileList   m_tileList;
     unsigned             m_liveTime;
     unsigned             m_prescaled;
     unsigned             m_discarded;
     sent_condArr         m_sent_condArr;
     unsigned             m_triggerTime;
     GemOverlayOnePpsTime m_onePpsTime;
     unsigned short       m_deltaEventTime;
     unsigned short       m_deltaWindOpenTime;
     unsigned short       m_missed;

 };


inline void GemOverlay::initTrigger(unsigned short tkr, unsigned short roi,
               unsigned short calLE, unsigned short calHE, 
               unsigned short cno, unsigned short conditionSummary, 
               unsigned short missed,
               const GemOverlayTileList &tileList) 
{
     m_tkrVector = tkr;
     m_roiVector = roi;
     m_cal_LE_Vector = calLE;
     m_cal_HE_Vector = calHE;
     m_cno_Vector = cno;
     m_conditionSummary = conditionSummary;
     m_missed = missed;
     m_tileList = tileList;
}

inline void GemOverlay::setTkrVector( unsigned short tkr)
{
    m_tkrVector = tkr;
}

/*
    inline void GemOverlay::initSummary(unsigned liveTime, unsigned prescaled, 
                  unsigned discarded, unsigned sent, unsigned triggerTime,
                  const GemOverlayOnePpsTime &time, unsigned short deltaEvtTime, unsigned short deltaWindOpenTime) {
        m_liveTime = liveTime;
        m_prescaled = prescaled;
        m_discarded = discarded;  
        m_sent_condArr.sent = sent;
        m_triggerTime = triggerTime;
        m_onePpsTime = time;
        m_deltaEventTime = deltaEvtTime;
        m_deltaWindOpenTime = deltaWindOpenTime;
    }
*/

inline void GemOverlay::initSummary(unsigned liveTime, unsigned prescaled, 
              unsigned discarded, GemOverlayDataCondArrivalTime condArr, unsigned triggerTime,
              const GemOverlayOnePpsTime &time, unsigned short deltaEvtTime, unsigned short deltaWindOpenTime) {
    m_liveTime = liveTime;
    m_prescaled = prescaled;
    m_discarded = discarded;  
    m_sent_condArr.condArr.init(condArr.condArr());
    m_triggerTime = triggerTime;
    m_onePpsTime = time;
    m_deltaEventTime = deltaEvtTime;
    m_deltaWindOpenTime = deltaWindOpenTime;
}

inline void GemOverlay::clear() {
    m_tkrVector = 0;
    m_roiVector = 0;  
    m_cal_LE_Vector = 0;
    m_cal_HE_Vector = 0;
    m_cno_Vector = 0;
    m_conditionSummary = 0;
    m_tileList.clear();
    m_liveTime = 0;
    m_prescaled = 0;
    m_discarded = 0;
    m_sent_condArr.sent = 0;
    m_sent_condArr.condArr.init(0);
    m_triggerTime = 0;
    m_onePpsTime.clear();
    m_deltaEventTime = 0;
    m_deltaWindOpenTime = 0;
    m_missed = 0;
}



inline std::ostream& GemOverlay::fillStream(std::ostream &s) const{
    s << "GEM:" <<std::endl;
    s << "ROI vector = 0x" << std::hex << std::setw(4) << std::setfill('0') << getRoiVector() << std::endl;
    s << "TKR vector = 0x" << std::hex << std::setw(4) << getTkrVector() << std::endl;
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
    return s;
}


inline std::ostream& operator<<(std::ostream &s, const GemOverlay &obj){
    return obj.fillStream(s);
}


/*!
* \brief  Return the condition arrival time in 50 ns ticks of the external trigger bit
* \return The condition arrival time in 50 ns ticks of the external trigger bit
*/
inline unsigned short  GemOverlayDataCondArrivalTime::external() const
{ 
  return (m_condArr >> 25) & ((1 << 5) -1);
}

/*!
* \brief  Return the condition arrival time in 50 ns ticks of the cno trigger bit
* \return The condition arrival time in 50 ns ticks of the cno trigger bit
*/
inline unsigned short  GemOverlayDataCondArrivalTime::cno()      const
{ 
  return (m_condArr >> 20) & ((1 << 5) -1);
}

/*!
* \brief  Return the condition arrival time in 50 ns ticks of the calHE trigger bit
* \return The condition arrival time in 50 ns ticks of the calHE trigger bit
*/
inline unsigned short  GemOverlayDataCondArrivalTime::calHE()    const
{ 
  return (m_condArr >> 15) & ((1 << 5) -1);
}

/*!
* \brief  Return the condition arrival time in 50 ns ticks of the calLE trigger bit
* \return The condition arrival time in 50 ns ticks of the calLE trigger bit
*/
inline unsigned short  GemOverlayDataCondArrivalTime::calLE()    const
{ 
  return (m_condArr >> 10) & ((1 << 5) -1);
}

/*!
* \brief  Return the condition arrival time in 50 ns ticks of the tkr trigger bit
* \return The condition arrival time in 50 ns ticks of the tkr trigger bit
*/
inline unsigned short  GemOverlayDataCondArrivalTime::tkr()      const
{ 
  return (m_condArr >>  5) & ((1 << 5) -1);
}

/*!
* \brief  Return the condition arrival time in 50 ns ticks of the roi trigger bit
* \return The condition arrival time in 50 ns ticks of the roi trigger bit
*/
inline unsigned short  GemOverlayDataCondArrivalTime::roi()      const
{ 
  return (m_condArr     ) & ((1 << 5) -1);
}



}//namespace GemOverlay


#endif
