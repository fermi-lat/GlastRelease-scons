
#ifndef IGcrSelectTool_h
#define IGcrSelectTool_h

#include "GaudiKernel/IAlgTool.h"
#include "Event/MonteCarlo/McParticle.h"

/**   
* @class IGcrSelectTool
*
* Base class for searching Cal data for GCRs
*
* 
*/

static const InterfaceID IID_IGcrSelectTool("IGcrSelectTool",1,0) ;

class IGcrSelectTool : virtual public IAlgTool 
{
  public:
    // retrieve Gaudi interface ID
    static const InterfaceID& interfaceID() { return IID_IGcrSelectTool; }

    //! main method
    
    virtual StatusCode selectGcrXtals()=0;
    
 } ;

#endif



