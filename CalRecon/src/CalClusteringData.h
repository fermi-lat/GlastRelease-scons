
#ifndef __CalClusteringData_H
#define __CalClusteringData_H 1

#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/IMessageSvc.h"
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "geometry/Point.h"

/**   
* @class CalClusteringData
*
* Data shared by all clustering actors, such as clustering algorithms
* and energy correction tools.
*
* $Header$
*/

class CalClusteringData
 {

  public:

    //! constructor
	CalClusteringData( ISvcLocator * ) ;

    //! update what depends on event
    void beginEvent() ;
    
    // generic services
    StatusCode getStatus() const
     { return m_status ; }
    IDataProviderSvc * getEventSvc() const
     { return m_eventSvc ; }
    IGlastDetSvc * getDetSvc() const
     { return m_detSvc ; }

    // cal data
    int getCalNLayers() const
     { return m_calNLayers ; }
    double getCalCsIWidth() const
     { return m_calCsIWidth ; }
    double getCalCsIHeight() const
     { return m_calCsIHeight ; }

    // tkr data
    int getTkrNVertices() const
     { return m_tkrNVertices ; }
    const Vector & getTkrFrontVertexDirection() const
     { return m_tkrFrontVertexDirection ; }
    const Point & getTkrFrontVertexPosition() const
     { return m_tkrFrontVertexPosition ; }

    // other
    void setSlope( double slope )
     { m_slope = slope ; }
    double getSlope() const
     { return m_slope ; }

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
    
    // tkr information
    int m_tkrNVertices ;
    Vector m_tkrFrontVertexDirection ;
    Point m_tkrFrontVertexPosition ;
    
    // other
    double m_slope;

 } ;

#endif



