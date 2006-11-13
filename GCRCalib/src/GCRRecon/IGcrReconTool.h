
#ifndef IGcrReconTool_h
#define IGcrReconTool_h

#include "GaudiKernel/IAlgTool.h"
#include "Event/MonteCarlo/McParticle.h"

/**   
* @class IGcrReconTool
*
* Base class for searching Cal data for GCRs
*
* 
*/

static const InterfaceID IID_IGcrReconTool("IGcrReconTool",1,0) ;

class IGcrReconTool : virtual public IAlgTool 
{
  public:
    // retrieve Gaudi interface ID
    static const InterfaceID& interfaceID() { return IID_IGcrReconTool; }

    //! main method
    
    virtual StatusCode findGcrXtals(bool useMcDir)=0;
    
 } ;

#endif



