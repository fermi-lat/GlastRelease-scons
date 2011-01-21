
#ifndef ICalClusterFiller_h
#define ICalClusterFiller_h

#include "Event/Recon/CalRecon/CalXtalRecData.h"
#include <CalRecon/ICalClusteringTool.h>

/**   
* @class ICalClusterFiller
*
* $Header$
*/


class ICalClusterFiller 
{
public:
    //* Defines the method for filling cluster info into CalCluster TDS objects
    virtual Event::CalCluster* fillClusterInfo(const XtalDataList* xtalVec) = 0;
};

#endif
        
        
        
