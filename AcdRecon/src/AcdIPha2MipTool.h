
#ifndef __AcdIPha2MipTool_H
#define __AcdIPha2MipTool_H 1

#include "GaudiKernel/IAlgTool.h"

#include "Event/Digi/AcdDigi.h"
#include "Event/Recon/AcdRecon/AcdHit.h"

#include "../AcdRecon/AcdReconStruct.h"

/**   
 * @class AcdIPha2MipTool
 *
 * @brief Gaudi interface for tool which converts raw PHA counts into MIP equivalent signals
 *
 * The actual conversions live in AcdUtil::AcdCalibFuncs
 * This code just loops over the objects and calls the relavent conversions routines
 * and makes new objects with the new calibrated quantities
 *
 * @author Eric Charles
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
  
  /// retrieve Gaudi interface ID
  static const InterfaceID& interfaceID()
    { return IID_AcdIPha2MipTool; }
  
  /// C'tor is trivial
  AcdIPha2MipTool() {}

  /// D'tor is trivial
  virtual ~AcdIPha2MipTool() {}
  
  
  /**
   * @brief Make collection of hits and fill the hit map
   * 
   * @param digiCol all the input AcdDigis
   * @param periodicEvent is the event a periodic trigger or not
   * @param hitCol collection to be filled with output AcdHit objects
   * @param hitMap hap to be filled with output AcdHit objects, maped by AcdId
   * @return Success or Failure
  **/
  virtual StatusCode makeAcdHits ( const Event::AcdDigiCol& digiCol,
				   bool periodicEvent, 
				   unsigned gemDeltaEventTime,
				   Event::AcdHitCol& hitCol,
				   AcdRecon::AcdHitMap& hitMap) = 0;
  
  /**
   * @brief Make a single hit
   * 
   * @param digi input AcdDigi
   * @param periodicEvent is the event a periodic trigger or not
   * @param hit newly made output AcdHit object
   * @return Success or Failure
  **/
  virtual StatusCode makeAcdHit ( const Event::AcdDigi& digi,
				  bool periodicEvent, 
				  Event::AcdHit*& hit) = 0;
 

} ;

#endif



