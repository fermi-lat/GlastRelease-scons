#ifndef AcdHeader_H
#define AcdHeader_H

#include "SubSystemHeader.h"

/*! \class AcdHeader
\brief Root class to contain the ACD header data for
one event. 
*/
class AcdHeader: public SubSystemHeader {

public:
    
    AcdHeader();

    virtual ~AcdHeader();

    virtual void Clean(Option_t *option="");


private:

    ClassDef(AcdHeader,1)
};

#endif

