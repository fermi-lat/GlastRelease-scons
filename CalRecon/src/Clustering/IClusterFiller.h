
#ifndef IClusterFiller_h
#define IClusterFiller_h

#include "Event/Recon/CalRecon/CalXtalRecData.h"
#include "ICalClusteringTool.h"

/**   
* @class IClusterFiller
*
* $Header$
*/


class IClusterFiller 
{
public:
    //* Defines the method for filling cluster info into CalCluster TDS objects
    virtual Event::CalCluster* fillClusterInfo(const XtalDataVec* xtalVec) = 0;
};

#endif
	
	
	
