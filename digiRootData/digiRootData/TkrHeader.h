#ifndef TkrHeader_H
#define TkrHeader_H

#include "SubSystemHeader.h"

/*! \class TkrHeader
\brief Root class to contain the TKR header data for
one event. 
*/
class TkrHeader: public SubSystemHeader {

public:
    
    TkrHeader();

    virtual ~TkrHeader();

    virtual void Clean(Option_t *option="");


private:

    ClassDef(TkrHeader,1)
};

#endif

