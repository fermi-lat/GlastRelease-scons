
#ifndef MomentsClusterInfo_h
#define MomentsClusterInfo_h

#include "Event/Recon/CalRecon/CalXtalRecData.h"
#include "Event/Recon/CalRecon/CalCluster.h"
#include "CalRecon/ICalClusterFiller.h"
#include "CalRecon/ICalReconSvc.h"
#include "CalMomentsAnalysis.h"

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


class MomentsClusterInfo : virtual public ICalClusterFiller
{
public:
    MomentsClusterInfo(const ICalReconSvc* calReconSvc, double transScaleFactor = 1.5);
    virtual ~MomentsClusterInfo() {};
    
    Event::CalCluster* fillClusterInfo(const XtalDataVec* xtalVec);
private:
    // Use this to fill the layer data, returns total energy
    double fillLayerData(const XtalDataVec* xtalVec, Event::CalCluster* cluster);

    // Use this to fill the moments information
    void fillMomentsData(const XtalDataVec* xtalVec, Event::CalCluster* cluster, double energy);

    //! package service
    const ICalReconSvc* m_calReconSvc;

    Point               m_p0;
    int                 m_calnLayers;
    CalMomentsDataVec   m_dataVec;

    /// Transverse RMS scale factor
    double m_transScaleFactor;
};

#endif
	
	
	
