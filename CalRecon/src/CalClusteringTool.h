
#ifndef __CalClusteringTool_H
#define __CalClusteringTool_H 1

#include "CalIClusteringTool.h"
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/ISvcLocator.h"
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "Event/Recon/CalRecon/CalXtalRecData.h"
#include "Event/Recon/CalRecon/CalCluster.h"

/**   
* @class CalClusteringTool
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


class CalClusteringTool : public CalIClusteringTool,  public AlgTool {
	
  public:
    
    CalClusteringTool
     ( const std::string & type, 
       const std::string & name,
       const IInterface * parent ) ;
    virtual ~CalClusteringTool() ;
    
	/// @brief Intialization of the tool
    virtual StatusCode initialize();

    /// @brief Default cluster finding framework
    virtual StatusCode findClusters(
        const Event::CalXtalRecCol *,
        Event::CalClusterCol *
      ) ;

  protected:
    
    //! useful types
    typedef  std::vector<Event::CalXtalRecData *> XtalDataVec ;
    typedef  std::vector<XtalDataVec *> XtalDataVecVec ;

	//! crystal width
    double m_CsIWidth;

    //! crystal height
    double m_CsIHeight;

    //! number of layers
    int m_CalnLayers;

    //! Collect CalXtalRecData pointers
    virtual void getXtals( const Event::CalXtalRecCol *, XtalDataVec & xtals ) ;

    //! Distinguish sets of related xtals
    virtual void makeSets( const XtalDataVec & xtals, XtalDataVecVec & clusters ) =0 ;

    //! Construct CalClusters from sets of related xtals
    virtual void setClusters( const XtalDataVecVec &, Event::CalClusterCol * ) ;

    // calculate the direction of the shower from the hits
    virtual Vector fitDirection
     ( std::vector<Vector> pos,
       std::vector<Vector> sigma2 ) ;
	
 } ;

#endif
	
	
	
