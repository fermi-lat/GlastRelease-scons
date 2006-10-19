#ifndef ENUMS_EVENTFLAGS_H
#define ENUMS_EVENTFLAGS_H

namespace enums {

    
    typedef enum {
            GOOD = 0,
            EVTSEQ = 1,
            TKRRECON = 2,
            PACKETERROR = 4,
            SUMMARYERROR = 8,
            TRGPARITYERROR = 16,
            BADLDFSTATUS = 32,
            GTRCPHASE = 64,
            GTFEPHASE = 128,
            GTCCFIFO = 256,
            GTCCHDRPARITY = 512,
            GTCCWCPARITY = 1024,
            GTRCSUMMARY = 2048,
            GTCCDATAPARITY = 4096,
            GTCCTIMEOUT = 8192,
            GCCCERROR = 16384,
            GTCCERROR = 32768,
            PHASEERROR = 65536,
            TIMEOUTERROR = 131072,
            MISSINGCALREADOUT= 262144,
            LDFWARNING = 524288
        } EventFlags;


}//namespace enums


#endif
