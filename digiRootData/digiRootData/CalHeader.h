#ifndef CalHeader_H
#define CalHeader_H

#include "SubSystemHeader.h"

/*! \class CalHeader
\brief Root class to contain the CAL header subpacket data for
one event. 
*/
class CalHeader: public SubSystemHeader {

public:
    
    CalHeader();

    virtual ~CalHeader();

    virtual void Clean(Option_t *option="");

private:
    
    ClassDef(CalHeader,1)
};

#endif

