
#ifndef __ICalReconSvc_H
#define __ICalReconSvc_H 1

#include "Event/Recon/CalRecon/CalXtalRecData.h"
#include "Event/Recon/CalRecon/CalCluster.h"
#include "Event/Recon/TkrRecon/TkrVertex.h"
#include "Event/Recon/TkrRecon/TkrTrack.h"
//#include "geometry/Point.h"
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/IInterface.h"

static const InterfaceID IID_ICalReconSvc("ICalReconSvc",1,0) ;

/**   
* @class ICalReconSvc
*
* Interface to data and features shared by all CalRecon actors.
*
* $Header$
*/

class ICalReconSvc : public virtual IInterface
 {

  public:

    //! retrieve Gaudi interface ID
    static const InterfaceID& interfaceID()
     { return IID_ICalReconSvc ; }
//    const IID & type() const
//     { return IID_ICalReconSvc ; }

   //! update what depends on event state
    virtual void reviewEvent() =0 ;
    
    // generic services
    virtual StatusCode getStatus() const =0 ;
    virtual IDataProviderSvc * getEventSvc() const =0 ;
    virtual IGlastDetSvc * getDetSvc() const =0 ;

    // cal parameters
    virtual int getCalNLayers() const =0 ;
    virtual double getCalCsIWidth() const =0 ;
    virtual double getCalCsIHeight() const =0 ;

    // cal event data
    virtual Event::CalXtalRecCol * getXtalRecs() =0 ;
    virtual Event::CalClusterCol * getClusters() =0 ;

    // tkr event data
    virtual int getTkrNVertices() const =0 ;
    virtual const Event::TkrVertex * getTkrFrontVertex() const =0 ;
    virtual int getTkrNTracks() const =0 ;
    virtual const Event::TkrTrack * getTkrFrontTrack() const =0 ;

    // utilities
    virtual double getSlope( Event::CalCluster * cluster ) const =0 ;

 } ;

#endif



