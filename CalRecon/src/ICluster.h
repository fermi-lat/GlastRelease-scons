
#ifndef __ICluster_H
#define __ICluster_H 1

#include "GaudiKernel/IAlgTool.h"
#include "Event/Recon/CalRecon/CalXtalRecData.h"
#include "Event/Recon/CalRecon/CalCluster.h"
#include "geometry/Vector.h"

/**   
* @class ICluster
*
* Base class for clustering tools
*
*
* $Header$
*/

static const InterfaceID IID_ICluster("ICluster", 1 , 0);

class ICluster : virtual public IAlgTool {

public:
    ICluster() {;};
    // retrieve interface ID
    static const InterfaceID& interfaceID() { return IID_ICluster; }
    //! constructor
    //! destructor
    virtual ~ICluster() {}; 

    virtual StatusCode initialize()=0;


    // worker function for finding clusters
    virtual StatusCode findClusters(Event::CalXtalRecCol* calXtalRecCol)=0;

    virtual StatusCode execute()=0;

    virtual StatusCode finalize()=0;

    virtual void setClusterCol(Event::CalClusterCol* calClusterCol)=0;


protected:
    virtual Vector Fit_Direction(std::vector<Vector> pos,
        std::vector<Vector> sigma2,
        int nlayers)=0;



};

#endif



