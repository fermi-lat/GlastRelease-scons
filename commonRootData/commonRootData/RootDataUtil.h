#ifndef ROOTDATAUTIL_H
#define ROOTDATAUTIL_H

#include "Exception.h"

namespace rootdatautil {
    extern void CheckBounds(long nIndex, long nLower, long nUpper);
    extern void CheckBufSize(unsigned long nBufSize, unsigned long nCorrectSize); 
};

#endif