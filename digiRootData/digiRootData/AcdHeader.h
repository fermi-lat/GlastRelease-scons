#ifndef AcdHeader_H
#define AcdHeader_H

#include "SubSystemHeader.h"

/*! \class AcdHeader
\brief ROOT class to contain the ACD header data for
one event. 
*/
class AcdHeader: public SubSystemHeader {

public:
    
    enum {
        L1T_Wait_Bit = 0x00100000,
        Cpu_Busy_Bit = 0x00080000,
        Fifo_Full_Trap_Bit = 0x00040000,
        Fifo_Full_Bit = 0x00020000,
        Fifo_Half_Trap_Bit = 0x00010000,
        Fifo_Half_Bit = 0x00008000,
        Acd_Busy_Bit =  0x00004000
    } DeadTimeCause;

    AcdHeader();

    virtual ~AcdHeader();

    virtual void Clean(Option_t *option="");

    bool isL1tWait() { return (m_deadtime_and_cause & L1T_Wait_Bit); };
    bool isCpuBusy() { return (m_deadtime_and_cause & Cpu_Busy_Bit); };
    bool isFifoFullTrap() { return (m_deadtime_and_cause & Fifo_Full_Trap_Bit); };
    bool isFifoFull() { return (m_deadtime_and_cause & Fifo_Full_Bit); };
    bool isFifoHalfTrap() { return (m_deadtime_and_cause & Fifo_Half_Trap_Bit); };
    bool isAcdBusy() { return (m_deadtime_and_cause & Acd_Busy_Bit); };


private:


    ClassDef(AcdHeader,1)
};

#endif

