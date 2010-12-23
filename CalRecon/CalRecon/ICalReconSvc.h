
#ifndef __ICalReconSvc_H
#define __ICalReconSvc_H 1

#include "Event/Recon/CalRecon/CalXtalRecData.h"
#include "Event/Recon/CalRecon/CalCluster.h"
#include "Event/Recon/TkrRecon/TkrVertex.h"
#include "Event/Recon/TkrRecon/TkrTrack.h"
#include "Event/Digi/CalDigi.h"
//#include "geometry/Point.h"
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/IInterface.h"
#include <string>

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
    static const InterfaceID& interfaceID() { return IID_ICalReconSvc ; }

    //! register errors
    virtual StatusCode handleError
      ( const std::string & catcherName,
        const std::string & comment ) =0 ;

    // generic services
    virtual IDataProviderSvc*     getEventSvc()     const = 0 ;
    virtual IGlastDetSvc*         getDetSvc()       const = 0 ;

    // cal parameters
    virtual int                   getCalNLayers()   const = 0 ;
    virtual double                getCalCsIWidth()  const = 0 ;
    virtual double                getCalCsIHeight() const = 0 ;
    virtual double                getCalCsILength() const = 0 ;
    virtual double                getCaltowerPitch() const = 0 ;
    virtual bool                  getCalFlightGeom() const = 0 ;

    // MomentsCLusterInfo parameters.
    virtual double getMciZeroSupprEnergy()           const = 0 ;
    virtual double getMciXtalsTruncFrac()            const = 0 ;
    virtual double getMciEneMomTruncFrac()           const = 0 ;

    // Quantities for the moments analysis.
    virtual double getMaTransScaleFactor()           const = 0 ;
    virtual double getMaTransScaleFactorBoost()      const = 0 ;
    virtual double getMaCoreRadius()                 const = 0 ;

    // cal event data
    virtual Event::CalXtalRecCol* getXtalRecs()           = 0 ;
    virtual Event::CalDigiCol *getDigis()                 = 0 ;
    virtual Event::CalClusterCol* getClusters()           = 0 ;
};

#endif



