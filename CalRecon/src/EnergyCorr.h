
#ifndef __EnergyCorr_H
#define __EnergyCorr_H 1

#include "IEnergyCorr.h"
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/ISvcLocator.h"
#include "Event/Recon/CalRecon/CalCluster.h"

/**   
* @class EnergyCorr
*
* Base class for energy corrections, containing common member data
*
* $Header$
*/

class EnergyCorr :  public IEnergyCorr, public AlgTool {
	
public:
    
    //! constructor
    EnergyCorr( const std::string & type, 
                const std::string & name, 
                const IInterface * parent)
     : AlgTool(type,name,parent)
     {}
     
    //! destructor
    virtual ~EnergyCorr()
     {} 
    
    double getEnergyCorr()
     { return m_energyCorr ; }

    void setEnergyCorr( double energyCorr )
     { m_energyCorr = energyCorr ; }
	
private:
	
    // energy correction
    double m_energyCorr ;
    
} ;

#endif
	
	
	
