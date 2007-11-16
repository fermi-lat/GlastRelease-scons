/** @file IInterleaveMap.h

    @brief Defines the IInterleave Map

$Header$

*/

#ifndef InterleaveMap_h
#define InterleaveMap_h

#include <string>
#include <vector>

/** @class InterleaveMap
    @brief Defines the structure of the "interleave map" which relates interleave events in the 
           output ntuple back to their original inpupt files
    @author Tracy Usher
*/

class InterleaveMap 
{
public:

    InterleaveMap() {}
    ~InterleaveMap() {}

public:
    // Below are the interleave map variables
    int    run;
    int    curEvent;   ///< current run, event
    int    irun;
    int    ievent; ///< interleaved run, event
    char   type[40];    ///< background type
    double value; ///< current value
};

class IInterleaveMap
{
public:
    virtual InterleaveMap* getInterleaveMap()     = 0;
    virtual void           saveInterleaveMapRow() = 0;
};


#endif
