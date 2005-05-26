
#ifndef MomentsClusterInfo_h
#define MomentsClusterInfo_h

#include "Event/Recon/CalRecon/CalXtalRecData.h"
#include "Event/Recon/CalRecon/CalCluster.h"
#include "IClusterFiller.h"
#include "Utilities/ICalReconSvc.h"

/**   
* @class MomentsClusterInfo
*
* Base class for clustering tools, containing member data and
* default code for the global algorithm, the preparation of
* a cluster from a set of crystals, and the computing of its
* direction.
* The only pure virtual method which needs to be implemented
* in a derived class is nextXtalsSet(), which is selecting the
* crystals to be grouped together.
*
* $Header$
*/


class MomentsClusterInfo : virtual public IClusterFiller
{
public:
    MomentsClusterInfo(const ICalReconSvc* calReconSvc);
   ~MomentsClusterInfo() {};
    
    Event::CalCluster* fillClusterInfo(const XtalDataVec* xtalVec);
private:
    // Use this to fill the layer data, returns total energy
    double fillLayerData(const XtalDataVec* xtalVec, Event::CalCluster* cluster);

    // Use this to fill the moments information
    void fillMomentsData(const XtalDataVec* xtalVec, Event::CalCluster* cluster, double energy);

    //! package service
    const ICalReconSvc* m_calReconSvc;

    const double        m_minFrac; 
    int                 m_calnLayers;
    const Point         m_p0;   // WBA: Spurious off-set for origin - ?????
};

#endif
	
	
	
