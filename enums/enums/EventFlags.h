#ifndef ENUMS_EVENTFLAGS_H
#define ENUMS_EVENTFLAGS_H

namespace enums {

    
    typedef enum {
            GOOD           =       0,
            EVTSEQ         =       1,
            TKRRECON       =       2,
            PACKETERROR    =       4,
            SUMMARYERROR   =       8,
            TRGPARITYERROR =    0x10, //16,
            BADLDFSTATUS   =    0x20, //32,
            GTRCPHASE      =    0x40, //64,
            GTFEPHASE      =    0x80, //128,
            GTCCFIFO       =   0x100, //256,
            GTCCHDRPARITY  =   0x200, //512,
            GTCCWCPARITY   =   0x400, //1024,
            GTRCSUMMARY    =   0x800, //2048,
            GTCCDATAPARITY =  0x1000, //4096,
            GTCCTIMEOUT    =  0x2000, //8192,
            GCCCERROR      =  0x4000, //16384,
            GTCCERROR      =  0x8000, //32768,
            PHASEERROR     = 0x10000, //65536,
            TIMEOUTERROR   = 0x20000, //131072,
            MISSINGCALREADOUT 
                           = 0x40000, //262144,
            LDFWARNING     = 0x80000  //524288
        } EventFlags;

namespace Gleam {
    typedef enum {
        GOOD = 0,
        TKRRECON = 1, 
        CALRECON = 2,
        ACDRECON = 4,
        ANALYSISNTUPLE = 8
     } GleamEventFlags; 

} //namespace Gleam

}//namespace enums


#endif
