
#ifndef __AcdIPha2MipTool_H
#define __AcdIPha2MipTool_H 1

#include "GaudiKernel/IAlgTool.h"

#include "Event/Digi/AcdDigi.h"
#include "Event/Recon/AcdRecon/AcdHit.h"

/**   
* @class AcdIPha2MipTool
*
* Base class for clustering tools
*
* $Header$
*/

static const InterfaceID IID_AcdIPha2MipTool("AcdIPha2MipTool",1,0) ;

class AcdIPha2MipTool : virtual public IAlgTool {

 public:
  
  // retrieve Gaudi interface ID
  static const InterfaceID& interfaceID()
    { return IID_AcdIPha2MipTool; }
  
  AcdIPha2MipTool() {}
  virtual ~AcdIPha2MipTool() {}
  
  //! main method
  virtual StatusCode makeAcdHits ( const Event::AcdDigiCol *,
				   Event::AcdHitCol * ) = 0;

} ;

#endif



