#ifndef LiveTime_H
#define LiveTime_H

#include "TObject.h"

/** @class LiveTime
* @brief data class to store live time packet data
* 
* The times in this class are actually 48 bit values.  There is no guarantee that 
* this the machine this class is run on supports 64 bit values...hence we define 
* overflow bits in a separate UShort_t.
*/


class LiveTime: public TObject {
private:
    /// number of triggers since start of the run
    UInt_t m_l1tCount;  

    /// 50 ns FPGA clock pulses without interruption to provide absolute uptime
    UInt_t m_l1tTime;
    UShort_t m_l1tTimeOverFlow;

    /// number of GPS receiver pulse per second pulses since power on
    UInt_t m_ppsCount;  

    /// value of L1tTime latched at the rising edge of the GPS receiver PPS signal
    UInt_t m_ppsTime;
    UShort_t m_ppsTimeOverFlow;
    
    /// counts 50 ns FPGA clock pulses while L1t enable is high AND L1T is high
    /// This indicates that a run is in progress and the instrument is ready for a trigger
    UInt_t m_liveTime;
    UShort_t m_liveTimeOverFlow;

    /// counts 50ns FPGA close pulses continously while L1t enable is high
    UInt_t m_runTime;
    UShort_t m_runTimeOverFlow;

    /// increments each time L1tCount Reset is given - at start of new run
    UInt_t m_runCount;

    UInt_t m_treqCapture;

    UInt_t m_ltcDeadWord;

public:

    LiveTime();
    virtual ~LiveTime();

    void Clean(Option_t *option="");
   
    void setL1TCount(UInt_t count) { m_l1tCount = count; };
    UInt_t getL1TCount() { return m_l1tCount; };

    void setL1TTime(ULong_t count) { m_l1tTime = count; };
    //! level one trigger time in counts
    ULong_t getL1TTime() { return m_l1tTime; };

    void setL1TTimeOverFlow(UShort_t over) { m_l1tTimeOverFlow = over; };
    //! returns the 16 bit overflow for the level 1 trigger time counter
    UShort_t getL1TTimeOverFlow() { return m_l1tTimeOverFlow; };

    void setPPSCount(UInt_t count) { m_ppsCount = count; };
    UInt_t getPPSCount() { return m_ppsCount; };

    void setPPSTime(ULong_t time) { m_ppsTime = time; };
    ULong_t getPPSTime() { return m_ppsTime; };

    void setPPSTimeOverFlow(UShort_t over) { m_ppsTimeOverFlow = over; };
    UShort_t getPPSTimeOverFlow() { return m_ppsTimeOverFlow; };

    void setLiveTime(ULong_t time) { m_liveTime = time; };
    //! returns the counts for the live time counter
    ULong_t getLiveTime() { return m_liveTime; };
     
    void setLiveTimeOverFlow(UShort_t over) { m_liveTimeOverFlow = over; };
    //! returns the overflow bits (16 bits) for live time counter
    UShort_t getLiveTimeOverFlow() { return m_liveTimeOverFlow; };

    void setRunTime(ULong_t time) { m_runTime = time; };
    //! returns the run time in counts
    ULong_t getRunTime() { return m_runTime; };

    void setRunTimeOverFlow(UShort_t over) { m_runTimeOverFlow = over; };
    //! returns the 16 bit overflow for the run time counter
    UShort_t getRunTimeOverFlow() { return m_runTimeOverFlow; };

    void setRunCount(UInt_t count) { m_runCount = count; };
    //! run number
    UInt_t getRunCount() { return m_runCount; };

    void setTREQCapture(UInt_t capture) { m_treqCapture = capture; };
    UInt_t getTREQCapture() { return m_treqCapture; };

    void setLTCDeadWord(UInt_t deadword) { m_ltcDeadWord = deadword; };
    UInt_t getLTCDeadWord() { return m_ltcDeadWord; };

    //ULong_t getLiveTimeFraction();


    
    ClassDef(LiveTime,1)       // Storage for live time data
}; 

#endif





