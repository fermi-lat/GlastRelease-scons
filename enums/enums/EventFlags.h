#ifndef ENUMS_EVENTFLAGS_H
#define ENUMS_EVENTFLAGS_H

namespace enums {

    
    typedef enum {
            GOOD = 0,
            EVTSEQ = 1,
            TKRRECON = 2,
            PACKETERROR = 4,
            SUMMARYERROR = 8,
            TRGPARITYERROR = 16
        } EventFlags;


}//namespace enums


#endif
