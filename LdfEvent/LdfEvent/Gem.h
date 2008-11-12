#ifndef LDFEVENT_LDFTIME_H
#define LDFEVENT_LDFTIME_H

#include <iostream>
#include <vector>

#include "GaudiKernel/Kernel.h"
#include "GaudiKernel/StreamBuffer.h"
#include "GaudiKernel/ContainedObject.h"
#include "GaudiKernel/ObjectVector.h"

#include "Event/TopLevel/Definitions.h"
#include "Event/TopLevel/EventModel.h"

/**
* @class Gem
* @brief TDS for GEM data
*/
static const CLID& CLID_LdfGem = InterfaceID("LdfGem", 1, 0);

namespace LdfEvent {

class GemTileList
{
public:
  GemTileList()  {}
  GemTileList (unsigned short xzm, unsigned short xzp, unsigned short yzm,
               unsigned short yzp, unsigned xy, unsigned short rbn, 
               unsigned short na) :
               m_xzp(xzp), m_xzm(xzm), m_yzp(yzp), m_yzm(yzm), m_xy(xy),
               m_rbn(rbn), m_na(na) { };

  GemTileList(const GemTileList &tileList) : 
               m_xzp(tileList.m_xzp), m_xzm(tileList.m_xzm), m_yzp(tileList.m_yzp), 
               m_yzm(tileList.m_yzm), m_xy(tileList.m_xy), m_rbn(tileList.m_rbn), 
               m_na(tileList.m_na) { };
               
  ~GemTileList() {}

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
		s << "XZM, XZP = 0x" << std::hex << std::setw(4) << xzm() << ", 0x"; 
        s << std::hex << xzp() << std::endl;
        s << "YZM, YZP = " << yzm() << ", " << yzp() << std::endl;
		s << "XY = " << std::setw(8) << xy() << std::endl;
		s << "RBN, NA = " << std::setw(4) << rbn() << ", " << na() << std::endl;
        return s;
    }

  unsigned short xzm()  const { return m_xzm; }
  unsigned short xzp()  const { return m_xzp; }
  unsigned short yzm()  const { return m_yzm; };
  unsigned short yzp()  const { return m_yzp; };
  unsigned       xy()   const { return m_xy; };
  unsigned short rbn()  const { return m_rbn; };
  unsigned short na()   const { return m_na; };

private:
  unsigned short m_xzp;
  unsigned short m_xzm;
  unsigned short m_yzp;
  unsigned short m_yzm;
  unsigned       m_xy;
  unsigned short m_rbn;
  unsigned short m_na;
};


/** @class GemOnePpsTime
 *  @brief This is a helper class for parsing the One PPS time structure.
 */
class GemOnePpsTime
{
public:
  GemOnePpsTime()  {}
  GemOnePpsTime(unsigned time, unsigned sec) { 
     m_timebase = time; 
     m_seconds = sec;
  };
  GemOnePpsTime(const GemOnePpsTime &time) : 
              m_timebase(time.m_timebase), m_seconds(time.m_seconds) {};

  ~GemOnePpsTime() {}

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

  unsigned timebase()  const { return m_timebase; };
  unsigned seconds()   const { return m_seconds; };

private:
  unsigned m_timebase;
  unsigned m_seconds;
};


class GemDataCondArrivalTime 
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
    class Gem : public DataObject{
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

        Gem() { clear(); };
        virtual ~Gem() { clear(); };
        virtual std::ostream& fillStream(std::ostream &s) const;
        friend std::ostream& operator << (std::ostream &s, const Gem &obj);


       // GEM doc v2.5 redefined sent() field as condition_arrival_times
        // We keep both accessors to this field for backward compatibility,
        // but overload the memory location with a union.
        union sent_condArr {
            unsigned sent;
            GemDataCondArrivalTime condArr;
        };

        /// Retrieve reference to class definition structure
        virtual const CLID& clID() const   { return Gem::classID(); }
        static const CLID& classID()       { return CLID_LdfGem; }

        void clear();

        void initTrigger(unsigned short tkr, unsigned short roi, 
                         unsigned short calLE,
                         unsigned short calHE, unsigned short cno, 
                         unsigned short conditionSummary, 
                         unsigned short missed,
                         const GemTileList &tileList);

        void setTkrVector( unsigned short tkr);

        void initSummary(unsigned liveTime, unsigned prescaled, 
                         unsigned discarded, GemDataCondArrivalTime condArr,
                         unsigned triggerTime, const GemOnePpsTime &time, 
                         unsigned short deltaEvtTime, unsigned short deltaWindOpenTime);
/*
        void initSummary(unsigned liveTime, unsigned prescaled, 
                         unsigned discarded, unsigned sent, 
                         unsigned triggerTime, const GemOnePpsTime &time, 
                         unsigned short deltaEvtTime, unsigned short deltaWindOpenTime); */


        unsigned short tkrVector() const { return m_tkrVector;};
        unsigned short roiVector() const { return m_roiVector;};
        unsigned short calLEvector() const { return m_cal_LE_Vector;};
        unsigned short calHEvector() const { return m_cal_HE_Vector;};
        unsigned short cnoVector() const { return m_cno_Vector;};
        unsigned short conditionSummary() const { return m_conditionSummary;};
        unsigned short missed() const { return m_missed; };
        const GemTileList& tileList() const { return m_tileList; };
        unsigned liveTime() const {return m_liveTime; };
        unsigned prescaled() const { return m_prescaled;};
        unsigned discarded() const { return m_discarded;};
        unsigned sent() const { return m_sent_condArr.sent; };
        GemDataCondArrivalTime condArrTime()   const { return m_sent_condArr.condArr; };
        unsigned triggerTime() const { return m_triggerTime; };
        const GemOnePpsTime& onePpsTime() const { return m_onePpsTime; };
        unsigned short deltaEventTime() const { return m_deltaEventTime; };
        unsigned short deltaWindowOpenTime() const { return m_deltaWindOpenTime; };

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

        unsigned short m_tkrVector;
        unsigned short m_roiVector;
        unsigned short m_cal_LE_Vector;
        unsigned short m_cal_HE_Vector;
        unsigned short m_cno_Vector;
        unsigned short m_conditionSummary;
        GemTileList m_tileList;
        unsigned  m_liveTime;
        unsigned  m_prescaled;
        unsigned  m_discarded;
        sent_condArr  m_sent_condArr;
        unsigned  m_triggerTime;
        GemOnePpsTime m_onePpsTime;
        unsigned short m_deltaEventTime;
        unsigned short m_deltaWindOpenTime;
        unsigned short m_missed;

    };


   inline void Gem::initTrigger(unsigned short tkr, unsigned short roi,
                  unsigned short calLE, unsigned short calHE, 
                  unsigned short cno, unsigned short conditionSummary, 
                  unsigned short missed,
                  const GemTileList &tileList) 
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

   inline void Gem::setTkrVector( unsigned short tkr)
   {
       m_tkrVector = tkr;
   }

/*
    inline void Gem::initSummary(unsigned liveTime, unsigned prescaled, 
                  unsigned discarded, unsigned sent, unsigned triggerTime,
                  const GemOnePpsTime &time, unsigned short deltaEvtTime, unsigned short deltaWindOpenTime) {
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

    inline void Gem::initSummary(unsigned liveTime, unsigned prescaled, 
                  unsigned discarded, GemDataCondArrivalTime condArr, unsigned triggerTime,
                  const GemOnePpsTime &time, unsigned short deltaEvtTime, unsigned short deltaWindOpenTime) {
        m_liveTime = liveTime;
        m_prescaled = prescaled;
        m_discarded = discarded;  
        m_sent_condArr.condArr.init(condArr.condArr());
        m_triggerTime = triggerTime;
        m_onePpsTime = time;
        m_deltaEventTime = deltaEvtTime;
        m_deltaWindOpenTime = deltaWindOpenTime;
    }

    inline void Gem::clear() {
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



    inline std::ostream& Gem::fillStream(std::ostream &s) const{
        s << "GEM:" <<std::endl;
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
        return s;
}


inline std::ostream& operator<<(std::ostream &s, const Gem &obj){
    return obj.fillStream(s);
}


/*!
* \brief  Return the condition arrival time in 50 ns ticks of the external trigger bit
* \return The condition arrival time in 50 ns ticks of the external trigger bit
*/
inline unsigned short  GemDataCondArrivalTime::external() const
{ 
  return (m_condArr >> 25) & ((1 << 5) -1);
}

/*!
* \brief  Return the condition arrival time in 50 ns ticks of the cno trigger bit
* \return The condition arrival time in 50 ns ticks of the cno trigger bit
*/
inline unsigned short  GemDataCondArrivalTime::cno()      const
{ 
  return (m_condArr >> 20) & ((1 << 5) -1);
}

/*!
* \brief  Return the condition arrival time in 50 ns ticks of the calHE trigger bit
* \return The condition arrival time in 50 ns ticks of the calHE trigger bit
*/
inline unsigned short  GemDataCondArrivalTime::calHE()    const
{ 
  return (m_condArr >> 15) & ((1 << 5) -1);
}

/*!
* \brief  Return the condition arrival time in 50 ns ticks of the calLE trigger bit
* \return The condition arrival time in 50 ns ticks of the calLE trigger bit
*/
inline unsigned short  GemDataCondArrivalTime::calLE()    const
{ 
  return (m_condArr >> 10) & ((1 << 5) -1);
}

/*!
* \brief  Return the condition arrival time in 50 ns ticks of the tkr trigger bit
* \return The condition arrival time in 50 ns ticks of the tkr trigger bit
*/
inline unsigned short  GemDataCondArrivalTime::tkr()      const
{ 
  return (m_condArr >>  5) & ((1 << 5) -1);
}

/*!
* \brief  Return the condition arrival time in 50 ns ticks of the roi trigger bit
* \return The condition arrival time in 50 ns ticks of the roi trigger bit
*/
inline unsigned short  GemDataCondArrivalTime::roi()      const
{ 
  return (m_condArr     ) & ((1 << 5) -1);
}



}//namespace LdfEvent


#endif
