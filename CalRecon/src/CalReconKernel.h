
#ifndef __CalReconKernel_H
#define __CalReconKernel_H 1

#include "Event/Recon/CalRecon/CalXtalRecData.h"
#include "Event/Recon/CalRecon/CalCluster.h"
#include "Event/Recon/TkrRecon/TkrVertex.h"
#include "Event/Recon/TkrRecon/TkrTrack.h"
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "geometry/Point.h"
#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/IMessageSvc.h"
#include "GaudiKernel/AlgTool.h"

static const InterfaceID IID_CalReconKernel("CalReconKernel",1,0) ;

/**   
* @class CalReconKernel
*
* Data shared by all CalRecon actors, such as clustering algorithms
* and energy correction tools.
*
* $Header$
*/

class CalReconKernel : public AlgTool
 {

  public:

    //! retrieve Gaudi interface ID
    static const InterfaceID& interfaceID()
     { return IID_CalReconKernel ; }

    //! constructor
	CalReconKernel
     ( const std::string & type, 
       const std::string & name,
       const IInterface * parent ) ;
    virtual ~CalReconKernel() ;

    //! init constant attributes
    virtual StatusCode initialize() ;
    
    //! update what depends on event state
    void reviewEvent() ;
    
    // generic services
    StatusCode getStatus() const
     { return m_status ; }
    IDataProviderSvc * getEventSvc() const
     { return m_eventSvc ; }
    IGlastDetSvc * getDetSvc() const
     { return m_detSvc ; }

    // cal parameters
    int getCalNLayers() const
     { return m_calNLayers ; }
    double getCalCsIWidth() const
     { return m_calCsIWidth ; }
    double getCalCsIHeight() const
     { return m_calCsIHeight ; }

    // cal event data
    Event::CalXtalRecCol * getXtalRecs()
     { return m_calXtalRecCol ; }
    Event::CalClusterCol * getClusters()
     { return m_calClusterCol ; }

    // tkr event data
    int getTkrNVertices() const
     { return m_tkrNVertices ; }
    const Event::TkrVertex * getTkrFrontVertex() const
     { return m_tkrFrontVertex ; }
    int getTkrNTracks() const
     { return m_tkrNTracks ; }
    const Event::TkrTrack * getTkrFrontTrack() const
     { return m_tkrFrontTrack ; }

    // utilities
    double getSlope( Event::CalCluster * cluster ) const
     { return (m_tkrSlope?m_tkrSlope:cluster->getDirection().z()) ; }

  private :
    
    // some pointers to services  
    StatusCode m_status ;
    IMessageSvc * m_messageSvc ;
    IDataProviderSvc * m_eventSvc ;
    IGlastDetSvc * m_detSvc; 
    
    //! CAL number of layers
    int m_calNLayers ;
    //! CAL crystal width
    double m_calCsIWidth ;
    //! CAL crystal height
    double m_calCsIHeight ;
    
    //! reconstructed data for crystals
    Event::CalXtalRecCol * m_calXtalRecCol ;
    //! the clusters list
    Event::CalClusterCol * m_calClusterCol ;

    // tkr information
    int m_tkrNVertices ;
    Event::TkrVertex * m_tkrFrontVertex ;
    double m_tkrSlope ;
    int m_tkrNTracks ;
    Event::TkrTrack * m_tkrFrontTrack ;

 } ;

#endif



