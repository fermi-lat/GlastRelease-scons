#ifndef ROOT_Gem_H
#define ROOT_Gem_H

#include "TObject.h"
#include "enums/GemConditionSummary.h"


class GemTileList : public TObject
{
public:
  GemTileList()  { Clear(); }
  GemTileList (UShort_t xzm, UShort_t xzp, UShort_t yzm,
               UShort_t yzp, UInt_t xy, UShort_t rbn, 
               UShort_t na) :
               m_xzp(xzp), m_xzm(xzm), m_yzp(yzp), m_yzm(yzm), m_xy(xy),
               m_rbn(rbn), m_na(na) { };

  GemTileList(const GemTileList &tileList) : TObject(tileList),
               m_xzp(tileList.m_xzp), m_xzm(tileList.m_xzm), m_yzp(tileList.m_yzp), 
               m_yzm(tileList.m_yzm), m_xy(tileList.m_xy), m_rbn(tileList.m_rbn), 
               m_na(tileList.m_na) { };
               
  ~GemTileList() {}

  void initialize (UShort_t xzm, UShort_t xzp, UShort_t yzm,
               UShort_t yzp, UInt_t xy, UShort_t rbn, 
               UShort_t na) { 
               m_xzp = xzp;
               m_xzm = xzm;
               m_yzp = yzp;
               m_yzm = yzm;
               m_xy = xy;
               m_rbn = rbn;
               m_na = na;
               };

  void Clear(Option_t *option="") {
               m_xzp = 0;
               m_xzm = 0;
               m_yzp = 0;
               m_yzm = 0;
               m_xy = 0;
               m_rbn = 0;
               m_na = 0;
  };

  void Print(Option_t *option="") const;

  /*@defgroup GemTileListGroup */
  /*@{*/
  UShort_t getXzm()  const { return m_xzm; }
  UShort_t getXzp()  const { return m_xzp; }
  UShort_t getYzm()  const { return m_yzm; };
  UShort_t getYzp()  const { return m_yzp; };
  UInt_t getXy()     const { return m_xy; };
  UShort_t getRbn()  const { return m_rbn; };
  UShort_t getNa()   const { return m_na; };
  /*@}*/

private:
  UShort_t m_xzp;
  UShort_t m_xzm;
  UShort_t m_yzp;
  UShort_t m_yzm;
  UInt_t m_xy;
  UShort_t m_rbn;
  UShort_t m_na;

  ClassDef(GemTileList,1) // Storage for GEM TileList
};


/** @class GemOnePpsTime
 *  @brief This is a helper class for parsing the One PPS time structure.
 */
class GemOnePpsTime : public TObject
{
public:
  GemOnePpsTime()  { Clear(); }
  GemOnePpsTime(UInt_t time, UInt_t sec) { 
     m_timebase = time; 
     m_seconds = sec;
  };
  GemOnePpsTime(const GemOnePpsTime &time) : TObject(time),
              m_timebase(time.m_timebase), m_seconds(time.m_seconds) {};

  ~GemOnePpsTime() {}

  void initialize(UInt_t timebase, UInt_t sec) { 
       m_timebase = timebase; 
       m_seconds = sec; 
  };

 
  void Print(Option_t *option="") const;

  void Clear(Option_t *option="") { 
     m_timebase = 0; 
     m_seconds = 0; 
  };

  /**@defgroup GemOnePpsTimeGroup GemOnePpsTime End-User Interface */
  /*@{*/
  UInt_t getTimebase()  const { return m_timebase; };
  UInt_t getSeconds()   const { return m_seconds; };
  /*@}*/

private:
  UInt_t m_timebase;
  UInt_t m_seconds;
  ClassDef(GemOnePpsTime,1) // Storage for GEM OnePPSTime
};

class GemCondArrivalTime : public TObject
{
public:

  GemCondArrivalTime() : TObject() { Clear(); };
  ~GemCondArrivalTime() { } ;

  void Clear() { m_condArr = 0; }

  void init(unsigned val) { m_condArr = val; };
/*!
* \brief  Return the condition arrival time in 50 ns ticks of the external trigg
er bit
* \return The condition arrival time in 50 ns ticks of the external trigger bit
*/
  UShort_t external()  const { return (m_condArr >> 25) & ((1 << 5) -1); }

/*!
* \brief  Return the condition arrival time in 50 ns ticks of the cno trigger bi
t
* \return The condition arrival time in 50 ns ticks of the cno trigger bit
*/
  UShort_t cno()       const { return (m_condArr >> 20) & ((1 << 5) -1); };


/*!
* \brief  Return the condition arrival time in 50 ns ticks of the calLE trigger
bit
* \return The condition arrival time in 50 ns ticks of the calLE trigger bit
*/
  UShort_t calLE()     const {  return (m_condArr >> 10) & ((1 << 5) -1); };


/*!
* \brief  Return the condition arrival time in 50 ns ticks of the calHE trigger
bit
* \return The condition arrival time in 50 ns ticks of the calHE trigger bit
*/
  UShort_t calHE()     const { return (m_condArr >> 15) & ((1 << 5) -1); }

/*!
* \brief  Return the condition arrival time in 50 ns ticks of the tkr trigger bi
t
* \return The condition arrival time in 50 ns ticks of the tkr trigger bit
*/
  UShort_t tkr()       const { return (m_condArr >>  5) & ((1 << 5) -1); }


/*!
* \brief  Return the condition arrival time in 50 ns ticks of the roi trigger bi
t
* \return The condition arrival time in 50 ns ticks of the roi trigger bit
*/
  UShort_t roi()       const { return (m_condArr     ) & ((1 << 5) -1); }


  /// Return fill packed word;
  UInt_t condArr() const { return m_condArr; };


private:
  UInt_t m_condArr;

  ClassDef(GemCondArrivalTime,1) // Storage for GEM Conditional Arrival Time 
};



/** @class Gem
 * @brief This stores GEM data.
 *
 * The Gem branch contains:
 * - TKR vector
 * - ROI vector 
 * - Condition Summary
 * - LiveTime
 * - prescaled
 * - Conditional Arrival Time
 *
 * For additional details concerning the GEM contribution, please refer to:
 * http://www-glast.slac.stanford.edu/IntegrationTest/ONLINE/docs/GEM.pdf
 * $Header$
 */

class Gem: public TObject {


public:

    Gem();
    Gem(const Gem& gem); 
    virtual ~Gem() { };

    void Clear(Option_t *option="");
 
    void Print(Option_t *option="") const;


    void initTrigger(UShort_t tkr, UShort_t roi, 
                     UShort_t calLE,
                     UShort_t calHE, UShort_t cno, 
                     UShort_t conditionSummary, 
                     const GemTileList &tileList);


    void initSummary(UInt_t liveTime, UInt_t prescaled, 
                     UInt_t discarded, UInt_t condArrTime, 
                     UInt_t triggerTime, const GemOnePpsTime &time, 
                     UShort_t deltaEvtTime, UShort_t deltaWindOpenTime);


	/** @ingroup GemGroup */
	/*@{*/
    UShort_t getTkrVector() const { return m_tkrVector;};
    UShort_t getRoiVector() const { return m_roiVector;};
    UShort_t getCalLeVector() const { return m_cal_Le_Vector;};
    UShort_t getCalHeVector() const { return m_cal_He_Vector;};
    UShort_t getCnoVector() const { return m_cno_Vector;};
    UShort_t getConditionSummary() const { return m_conditionSummary;};
    UShort_t getMissed() const { return (_conditionSummary >> 8) & (( 1 << 8 ) - 1); };
    const GemTileList& getTileList() const { return m_tileList; };
    UInt_t getLiveTime() const {return m_liveTime; };
    UInt_t getPrescaled() const { return m_prescaled;};
    UInt_t getDiscarded() const { return m_discarded;};
//    UInt_t getSent() const { return m_sent; };
    const GemCondArrivalTime& getCondArrTime()   const { return m_condArr; };
    UInt_t getTriggerTime() const { return m_triggerTime; };
    const GemOnePpsTime& getOnePpsTime() const { return m_onePpsTime; };
    UShort_t getDeltaEventTime() const { return m_deltaEventTime; };
    UShort_t getDeltaWindowOpenTime() const { return m_deltaWindOpenTime; };


    /// Method to query ROI bit in the condition summary word
    Bool_t getRoiSet() const { return( (m_conditionSummary & enums::ROI) != 0); };
	/// Method to query TKR bit in condition summary word
    Bool_t getTkrSet() const { return( (m_conditionSummary & enums::TKR) != 0); };
	/// Method to query Cal LE bit in condition summary word
    Bool_t getCalLeSet() const { return( (m_conditionSummary & enums::CALLE) != 0); };
	/// Method Cal HE bit in condition summary word
    Bool_t getCalHeSet() const { return( (m_conditionSummary & enums::CALHE) != 0); };
	/// Method to query CNO bit in condition summary word
    Bool_t getCnoSet() const { return( (m_conditionSummary & enums::CNO) != 0); };
	/// Method to query periodic bit in condition summary word
    Bool_t getPeriodicSet() const { return( (m_conditionSummary & enums::PERIODIC) != 0); };
	/// Query Solicited bit in condition summary word
    Bool_t getSolicitedSet() const { return( (m_conditionSummary & enums::SOLICITED) != 0); };
	/// Query External bit in condition summary word
    Bool_t getExternalSet() const { return( (m_conditionSummary & enums::EXTERNAL) != 0); };
	/*@}*/

    private:

    UShort_t m_tkrVector;
    UShort_t m_roiVector;
    UShort_t m_cal_Le_Vector;
    UShort_t m_cal_He_Vector;
    UShort_t m_cno_Vector;
    UShort_t m_conditionSummary;
    GemTileList m_tileList;
    UInt_t  m_liveTime;
    UInt_t m_prescaled;
    UInt_t  m_discarded;
    UInt_t  m_sent;
    UInt_t m_triggerTime;
    GemOnePpsTime m_onePpsTime;
    UInt_t  m_deltaEventTime;
    GemCondArrivalTime m_condArr;
    UShort_t m_deltaWindOpenTime;

    ClassDef(Gem,2) // Storage for GEM 
}; 
 
#endif


