
#ifndef __CalEnergyCorr_H
#define __CalEnergyCorr_H 1

#include "ICalEnergyCorr.h"
#include "ICalReconSvc.h"
#include "Event/Recon/CalRecon/CalCluster.h"
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/ISvcLocator.h"

/**   
* @class CalEnergyCorr
*
* Base class for energy corrections
*
* $Header$
*/

class CalEnergyCorr :  public ICalEnergyCorr, public AlgTool {
	
public:
    
    //! constructor
    CalEnergyCorr( const std::string & type, 
                const std::string & name, 
                const IInterface * parent)
     : AlgTool(type,name,parent), m_calReconSvc(0)
     {}
     
    virtual StatusCode initialize() ;

    //! calculate corrections on all cal clusters
    virtual StatusCode doEnergyCorr() ;
    
    //! destructor
    virtual ~CalEnergyCorr()
     {} 

  protected :
  
    //! calculate corrections on a given cal cluster
    virtual StatusCode doEnergyCorr( Event::CalCluster * ) =0 ;
    
    //! package service
    ICalReconSvc * m_calReconSvc ;
    
} ;

#endif
	
	
	
