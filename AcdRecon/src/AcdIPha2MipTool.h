
#ifndef __AcdIPha2MipTool_H
#define __AcdIPha2MipTool_H 1

#include "GaudiKernel/IAlgTool.h"

#include "Event/Digi/AcdDigi.h"
#include "Event/Recon/AcdRecon/AcdHit.h"

#include "../AcdRecon/AcdReconStruct.h"

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

  enum { AcceptMapBit_AMask = 0x1,
	 AcceptMapBit_BMask = 0x2,
	 VetoBit_AMask = 0x4,
	 VetoBit_BMask = 0x8,
	 CNO_AMask = 0x10,
	 CNO_BMask = 0x20 } MASKS;

public:
  
  // retrieve Gaudi interface ID
  static const InterfaceID& interfaceID()
    { return IID_AcdIPha2MipTool; }
  
  AcdIPha2MipTool() {}
  virtual ~AcdIPha2MipTool() {}
  
  
  /// @brief Make collection of hits and fill the hit map
  virtual StatusCode makeAcdHits ( const Event::AcdDigiCol&,
				   Event::AcdHitCol&,
				   AcdRecon::AcdHitMap&) = 0;
  
  /// @brief Make a single hit
  virtual StatusCode makeAcdHit ( const Event::AcdDigi&,
				  Event::AcdHit*& ) = 0;

} ;

#endif



