#ifndef CalHeader_H
#define CalHeader_H

#include "SubSystemHeader.h"

/*! \class CalHeader
\brief ROOT class to contain the CAL header subpacket data for
one event. 
*/
class CalHeader: public SubSystemHeader {

public:

   enum {
        L1T_Wait_Bit = 0x00100000,
        Cpu_Busy_Bit = 0x00080000,
        Fifo_Full_Trap_Bit = 0x00040000,
        Fifo_Full_Bit = 0x00020000,
        Fifo_Half_Trap_Bit = 0x00010000,
        Fifo_Half_Bit = 0x00008000,
        Cal_Busy_Bit =  0x00004000
   } DeadTimeCause;


    CalHeader();

    virtual ~CalHeader();

    virtual void Clean(Option_t *option="");

    bool isL1tWait() { return (m_deadtime_and_cause & L1T_Wait_Bit); };
    bool isCpuBusy() { return (m_deadtime_and_cause & Cpu_Busy_Bit); };
    bool isFifoFullTrap() { return (m_deadtime_and_cause & Fifo_Full_Trap_Bit); };
    bool isFifoFull() { return (m_deadtime_and_cause & Fifo_Full_Bit); };
    bool isFifoHalfTrap() { return (m_deadtime_and_cause & Fifo_Half_Trap_Bit); };
    bool isCalBusy() { return (m_deadtime_and_cause & Cal_Busy_Bit); };


private:

    
    ClassDef(CalHeader,1)
};

#endif

