
#ifndef __ClusteringTool_H
#define __ClusteringTool_H 1

#include "IClusteringTool.h"
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/ISvcLocator.h"
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "Event/Recon/CalRecon/CalXtalRecData.h"
#include "Event/Recon/CalRecon/CalCluster.h"

/**   
* @class ClusteringTool
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


class ClusteringTool : public IClusteringTool,  public AlgTool {
	
  public:
    
    ClusteringTool
     ( const std::string & type, 
       const std::string & name,
       const IInterface * parent ) ;
    virtual ~ClusteringTool() ;
    
	/// @brief Intialization of the tool
    virtual StatusCode initialize();

    /// @brief Default cluster finding framework
    virtual StatusCode findClusters(
        Event::CalXtalRecCol *,
        Event::CalClusterCol *
      ) ;

  protected:
    
    //! useful type
    typedef  std::vector<Event::CalXtalRecData*> xTalDataVec ;

	//! crystal width
    double m_CsIWidth;

    //! crystal height
    double m_CsIHeight;

    //! number of layers
    int m_CalnLayers;

    //! This should find the next highest energy cluster from CalXtalRecData pointers
    virtual xTalDataVec nextXtalsSet( xTalDataVec & xTalVec ) =0 ;

    //! This makes a CalCluster out of CalXtalRecData pointers
    virtual void makeCluster( xTalDataVec &, Event::CalClusterCol * ) ;

    // calculate the direction of the shower from the hits
    virtual Vector fitDirection
     ( std::vector<Vector> pos,
       std::vector<Vector> sigma2 ) ;
	
 } ;

#endif
	
	
	
