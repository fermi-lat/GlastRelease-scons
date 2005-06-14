
#ifndef IMipFindingTool_h
#define IMipFindingTool_h

#include "GaudiKernel/IAlgTool.h"

/**   
* @class IMipFindingTool
*
* Base class for searching Cal data for MIPs
*
* $Header$
*/

static const InterfaceID IID_IMipFindingTool("IMipFindingTool",1,0) ;

class IMipFindingTool : virtual public IAlgTool 
{
  public:
    // retrieve Gaudi interface ID
    static const InterfaceID& interfaceID() { return IID_IMipFindingTool; }

    //! main method
    virtual StatusCode findMIPCandidates() =0 ;
 } ;

#endif



