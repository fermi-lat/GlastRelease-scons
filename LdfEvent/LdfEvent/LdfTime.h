#ifndef LDFCONVERTER_LDFTIME_H
#define LDFCONVERTER_LDFTIME_H

#include <iostream>
#include <vector>

#include "GaudiKernel/Kernel.h"
#include "GaudiKernel/StreamBuffer.h"
#include "GaudiKernel/ContainedObject.h"
#include "GaudiKernel/ObjectVector.h"

#include "Event/TopLevel/Definitions.h"
#include "Event/TopLevel/EventModel.h"

/**
* @class EbfTime
* @brief TDS for time
*/
extern const CLID& CLID_LdfTime;
namespace LdfEvent {
    class LdfTime : public DataObject{
    public:
        virtual ~LdfTime();
        virtual std::ostream& fillStream(std::ostream &s) const;
        friend std::ostream& operator << (std::ostream &s, const LdfTime &obj);

        LdfTime();
        LdfTime(unsigned timeSec, unsigned timeNanoSec, unsigned upperPpcTimeBase, unsigned lowerPpcTimeBase);
        void initialize(unsigned timeSec, unsigned timeNanoSec, unsigned upperPpcTimeBase, unsigned lowerPpcTimeBase);

        unsigned timeSec() const { return m_timeSec; }
        unsigned timeNanoSec() const {return m_timeNanoSec; }
        unsigned upperPpcTimeBaseWord() const { return m_upperPpcTimeBase; };
        unsigned lowerPpcTimeBaseWord() const { return m_lowerPpcTimeBase; };

    private:
        unsigned m_timeSec;
        unsigned m_timeNanoSec;
        
        /// These two members form a 64 bit word that represents a "clock cycle number" from
        /// the PPC registers - this value divided by 16MHz is approximately the time since 
        /// the VxWorks box was turned on
        unsigned m_upperPpcTimeBase;
        unsigned m_lowerPpcTimeBase;

    };

    // the put methods, here for now
    inline LdfTime::LdfTime(){
        m_timeSec = 0;
        m_timeNanoSec = 0;
    }

    inline LdfTime::LdfTime(unsigned timeSec, unsigned timeNanoSec, 
        unsigned upperPpcTimeBase, unsigned lowerPpcTimeBase) {
        m_timeSec = timeSec;
        m_timeNanoSec = timeNanoSec;
        m_upperPpcTimeBase = upperPpcTimeBase;
        m_lowerPpcTimeBase = lowerPpcTimeBase;
    }

    inline void LdfTime::initialize(unsigned timeSec, unsigned timeNanoSec,
        unsigned upperPpcTimeBase, unsigned lowerPpcTimeBase) {
        m_timeSec = timeSec;
        m_timeNanoSec = timeNanoSec;
        m_upperPpcTimeBase = upperPpcTimeBase;
        m_lowerPpcTimeBase = lowerPpcTimeBase;
    }

    inline LdfTime::~LdfTime(){
    }


    inline std::ostream& LdfTime::fillStream(std::ostream &s) const{
        s<<"Time: " << m_timeSec << " Seconds " << m_timeNanoSec << " NanoSeconds " <<std::endl;
        return s;
    }


    inline std::ostream& operator<<(std::ostream &s, const LdfTime &obj){
        return obj.fillStream(s);
    }

}//namespace LdfEvent


#endif
