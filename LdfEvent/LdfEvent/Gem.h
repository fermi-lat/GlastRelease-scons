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
extern const CLID& CLID_LdfGem;
namespace LdfEvent {

class GemDataTileList
{
public:
  GemDataTileList()  {}
  GemDataTileList (unsigned short xzm, unsigned short xzp, unsigned short yzm,
               unsigned short yzp, unsigned xy, unsigned short rbn, 
               unsigned short na) :
               m_xzp(xzp), m_xzm(xzm), m_yzp(yzp), m_yzm(yzm), m_xy(xy),
               m_rbn(rbn), m_na(na) { };

  GemDataTileList(const GemDataTileList &tileList) : 
               m_xzp(tileList.m_xzp), m_xzm(tileList.m_xzm), m_yzp(tileList.m_yzp), 
               m_yzm(tileList.m_yzm), m_xy(tileList.m_xy), m_rbn(tileList.m_rbn), 
               m_na(tileList.m_na) { };
               
  ~GemDataTileList() {}

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

    inline std::ostream& Gem::fillStream(std::ostream &s) const{
        s << "tile list: " << std::endl;
        s << "XZM, XZP = 0x"" << std::hex << xzm() << ", 0x"; 
        s std::hex << xzy() << std::endl;
        s << "YZM, YZP = " << yzm() << ", " << yzp() << std::endl;
        s << "XY = " << xy() << std::endl;
        s << "RBN, NA = " << rbn() << ", " << na() << std::endl;
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


/** @class GemDataOnePpsTime
 *  @brief This is a helper class for parsing the One PPS time structure.
 */
class GemDataOnePpsTime
{
public:
  GemDataOnePpsTime()  {}
  GemDataOnePpsTime(unsigned time, unsigned sec) { 
     m_timebase = time; 
     m_seconds = sec;
  };
  GemDataOnePpsTime(const GemDataOnePpsTime &time) : 
              m_timebase(time.m_timebase), m_seconds(time.m_seconds) {};

  ~GemDataOnePpsTime() {}

  void init(unsigned timebase, unsigned sec) { 
       m_timebase = timebase; 
       m_seconds = sec; 
  };

  void clear() { 
     m_timebase = 0; 
     m_seconds = 0; 
  };

  inline std::ostream& Gem::fillStream(std::ostream &s) const{
      s << "OnePPSTime:" << endl;
      s << "seconds, timebase = " << m_seconds << ", " << m_timebase <<endl;
      return s;
  }

  unsigned timebase()  const { return m_timebase; };
  unsigned seconds()   const { return m_seconds; };

private:
  unsigned m_timebase;
  unsigned m_seconds;
};




    /** @class Gem
      * @brief Local storage of GEM data
      * $Header$
    */
    class Gem : public DataObject{
    public:
        Gem() { clear(); };
        virtual ~Gem() { clear(); };
        virtual std::ostream& fillStream(std::ostream &s) const;
        friend std::ostream& operator << (std::ostream &s, const Gem &obj);

        void clear();

        void initTrigger(unsigned short tkr, unsigned short roi, 
                         unsigned short calLE,
                         unsigned short calHE, unsigned short cno, 
                         unsigned short conditionSummary, 
                         const GemDataTileList &tileList);

        void initSummary(unsigned liveTime, unsigned prescaled, 
                         unsigned discarded, unsigned sent, 
                         unsigned triggerTime, const GemDataOnePpsTime &time, 
                         unsigned deltaEvtTime);


        unsigned short tkrVector() const { return m_tkrVector;};
        unsigned short roiVector() const { return m_roiVector;};
        unsigned short calLEvector() const { return m_cal_LE_Vector;};
        unsigned short calHEvector() const { return m_cal_HE_Vector;};
        unsigned short cnoVector() const { return m_cno_Vector;};
        unsigned short conditionSummary() const { return m_conditionSummary;};
        const GemDataTileList& tileList() const { return m_tileList; };
        unsigned liveTime() const {return m_liveTime; };
        unsigned prescaled() const { return m_prescaled;};
        unsigned discarded() const { return m_discarded;};
        unsigned sent() const { return m_sent; };
        unsigned triggerTime() const { return m_triggerTime; };
        const GemDataOnePpsTime& onePpsTime() const { return m_onePpsTime; };
        unsigned deltaEventTime() const { return m_deltaEventTime; };


    private:

        unsigned short m_tkrVector;
        unsigned short m_roiVector;
        unsigned short m_cal_LE_Vector;
        unsigned short m_cal_HE_Vector;
        unsigned short m_cno_Vector;
        unsigned short m_conditionSummary;
        GemDataTileList m_tileList;
        unsigned  m_liveTime;
        unsigned  m_prescaled;
        unsigned  m_discarded;
        unsigned  m_sent;
        unsigned  m_triggerTime;
        GemDataOnePpsTime m_onePpsTime;
        unsigned  m_deltaEventTime;

    };


   void Gem::initTrigger(unsigned short tkr, unsigned short roi,
                  unsigned short calLE, unsigned short calHE, 
                  unsigned short cno, unsigned short conditionSummary, 
                  const GemDataTileList &tileList) {
        m_tkrVector = tkr;
        m_roiVector = roi;
        m_cal_LE_Vector = calLE;
        m_cal_HE_Vector = calHE;
        m_cno_Vector = cno;
        m_conditionSummary = conditionSummary;
        m_tileList = tileList;
    }

    void Gem::initSummary(unsigned liveTime, unsigned prescaled, 
                  unsigned discarded, unsigned sent, unsigned triggerTime,
                  const GemDataOnePpsTime &time, unsigned deltaEvtTime) {
        m_liveTime = liveTime;
        m_prescaled = prescaled;
        m_discarded = discarded;  
        m_sent = sent;
        m_triggerTime = triggerTime;
        m_onePpsTime = time;
        m_deltaEventTime = deltaEvtTime;
    }

    void Gem::clear() {
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
        m_sent = 0;
        m_triggerTime = 0;
        m_onePpsTime.clear();
        m_deltaEventTime = 0;
    }



    inline std::ostream& Gem::fillStream(std::ostream &s) const{
        s << "GEM:" <<std::endl;
        s << "ROI vector = 0x" << std::hex << roiVector() << std::endl;
        s << "TKR vector = 0x" << std::hex << tkrVector() << std::endl;
        s << "CAL HE vector = 0x" << std::hex << m_cal_HE_Vector << std::endl;
        s << "CAL LE vector = 0x" << std::hex << m_cal_LE_Vector << std::endl;
        s << "Condition Summary = 0x" << std::hex << m_conditionSummary << std::endl;
        s << "CNO vector        = 0x" << std::hex << m_cno_Vector << std::endl;
        m_tileList.fillStream(s);
        s << "Live time         = 0x" << std::hex << m_liveTime << " = "; 
        s << m_liveTime << std::endl;
        s << "Prescaled         = 0x" << std::hex << m_prescaled << " = "; 
        s << m_prescaled << std::endl;
        s << "Discarded         = 0x" << std::hex << m_discarded << " = ";
        s << m_discarded << std::endl;
        s << "Sent              = 0x" << std::hex << m_sent << " = " ;
        s << m_sent << std::endl;
        m_onePpsTime.fillStream(s);
        s << "Delta event time = 0x" << std::hex << m_deltaEventTime << " = ";
        s << m_deltaEventTime << std::endl;
        return s;
    }


    inline std::ostream& operator<<(std::ostream &s, const Gem &obj){
        return obj.fillStream(s);
    }

}//namespace LdfEvent


#endif
