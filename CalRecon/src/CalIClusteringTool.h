
#ifndef __CalIClusteringTool_H
#define __CalIClusteringTool_H 1

#include "GaudiKernel/IAlgTool.h"
#include "Event/Recon/CalRecon/CalXtalRecData.h"
#include "Event/Recon/CalRecon/CalCluster.h"
#include "geometry/Vector.h"

/**   
* @class CalIClusteringTool
*
* Base class for clustering tools
*
* $Header$
*/

static const InterfaceID IID_CalIClusteringTool("CalIClusteringTool",1,0) ;

class CalIClusteringTool : virtual public IAlgTool {

  public:

    // retrieve Gaudi interface ID
    static const InterfaceID& interfaceID()
     { return IID_CalIClusteringTool; }

    CalIClusteringTool() {}
    virtual ~CalIClusteringTool() {}

    //! main method
    virtual StatusCode findClusters(
        Event::CalXtalRecCol *,
        Event::CalClusterCol *
      ) =0 ;

 } ;

#endif



