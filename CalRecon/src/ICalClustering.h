
#ifndef __ICalClustering_H
#define __ICalClustering_H 1

#include "ICalReconSvc.h"
#include "geometry/Vector.h"
#include "Event/Recon/CalRecon/CalXtalRecData.h"
#include "Event/Recon/CalRecon/CalCluster.h"
#include "GaudiKernel/IAlgTool.h"

/**   
* @class ICalClustering
*
* Base class for clustering tools
*
* $Header$
*/

static const InterfaceID IID_ICalClustering("ICalClustering",1,0) ;

class ICalClustering : virtual public IAlgTool {

  public:

    // retrieve Gaudi interface ID
    static const InterfaceID& interfaceID()
     { return IID_ICalClustering; }

    ICalClustering() {}
    virtual ~ICalClustering() {}

    //! main method
    virtual StatusCode findClusters() =0 ;

 } ;

#endif



