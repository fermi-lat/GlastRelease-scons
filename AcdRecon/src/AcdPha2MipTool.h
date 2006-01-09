
#ifndef __AcdPha2MipTool_H
#define __AcdPha2MipTool_H 1

#include "AcdIPha2MipTool.h"
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/ISvcLocator.h"


class MsgStream;
class HepMatrix;

/**   
* @class AcdPha2MipTool
*
* Class for extrapolating the Tracks to their intersections w/ the ACD
*
* $ $
*/


class AcdPha2MipTool : public AcdIPha2MipTool,  public AlgTool {
	
public:
    
  AcdPha2MipTool
  ( const std::string & type, 
    const std::string & name,
    const IInterface * parent ) ;
  virtual ~AcdPha2MipTool() ;
  
  /// @brief Intialization of the tool
  virtual StatusCode initialize();

  /// @brief Default cluster finding framework
  virtual StatusCode makeAcdHits ( const Event::AcdDigiCol *,
				   Event::AcdHitCol* );
  
protected:
  
  bool getCalibratedValues(const Event::AcdDigi& digi, float& mipsPmtA, float& mipsPmtB) const;

private:

  Event::AcdHitCol* output;

} ;

#endif
	
	
	
