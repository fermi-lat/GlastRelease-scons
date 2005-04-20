
#ifndef __CalEnergyCorr_H
#define __CalEnergyCorr_H 1

#include "ICalEnergyCorr.h"
#include "CalReconActor.h"
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

class CalEnergyCorr :  public ICalEnergyCorr, public AlgTool, protected CalReconActor {
	
public:
    
    //! constructor
    CalEnergyCorr( const std::string & type, 
                const std::string & name, 
                const IInterface * parent)
     : AlgTool(type,name,parent)
     {}
     
    virtual StatusCode initialize() ;

    //! calculate corrections on all cal clusters
    virtual StatusCode doEnergyCorr() ;
    
    //! calculate corrections on a given cal cluster
    virtual StatusCode doEnergyCorr( Event::CalCluster * ) =0 ;
    
    //! destructor
    virtual ~CalEnergyCorr()
     {} 
        
} ;

#endif
	
	
	
