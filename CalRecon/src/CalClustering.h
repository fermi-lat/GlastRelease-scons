
#ifndef __CalClustering_H
#define __CalClustering_H 1

#include "ICalClustering.h"
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "Event/Recon/CalRecon/CalXtalRecData.h"
#include "Event/Recon/CalRecon/CalCluster.h"
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/ISvcLocator.h"

/**   
* @class CalClustering
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


class CalClustering : public ICalClustering,  public AlgTool {
	
  public:
    
    CalClustering
     ( const std::string & type, 
       const std::string & name,
       const IInterface * parent ) ;
    virtual ~CalClustering() ;
    
	/// @brief Intialization of the tool
    virtual StatusCode initialize() ;

    /// @brief Default cluster finding framework
    virtual StatusCode findClusters() ;

  protected:
    
    //! package service
    ICalReconSvc * m_calReconSvc ;
    
    //! useful types
    typedef  std::vector<Event::CalXtalRecData *> XtalDataVec ;
    typedef  std::vector<XtalDataVec *> XtalDataVecVec ;

    //! Collect CalXtalRecData pointers
    virtual void getXtals( XtalDataVec & xtals ) ;

    //! THE METHOD TO IMPLEMENT: aggregate the sets of related xtals
    virtual void makeSets( const XtalDataVec & xtals, XtalDataVecVec & clusters ) =0 ;

    //! Construct CalClusters from sets of related xtals
    virtual void setClusters( const XtalDataVecVec & ) ;

    // calculate the direction of the shower from the hits
    virtual Vector fitDirection
     ( std::vector<Vector> pos,
       std::vector<Vector> sigma2 ) ;

 } ;

#endif
	
	
	
