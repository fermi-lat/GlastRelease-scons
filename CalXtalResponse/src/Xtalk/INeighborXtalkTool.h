#ifndef INeighborXtalkTool_h
#define INeighborXtalkTool_h
//  $Header$

/** @file
    @author Z.Fewtrell
*/


// LOCAL INCLUDES

// GLAST INCLUDES
#include "CalUtil/CalDefs.h"
#include "Event/Digi/CalDigi.h"


// EXTLIB INCLUDES
#include "GaudiKernel/IAlgTool.h"

// STD INCLUDES

static const InterfaceID IID_INeighborXtalkTool("INeighborXtalkTool", 1, 2);


/** 
    \author Z.Fewtrell

    \brief Abstract interface class for calculation of electronic crosstalk
    in GLAST Cal originating from neighboring CsI crystal detector element electronics.

    Generally speaking, cross-talk @ any point will depend on a static crosstalk coefficient table & a per-event signal map for the whole cal.
    Tool implementation will then model crosstalk measured @ any particular channel based on these tables.
*/
class INeighborXtalkTool : virtual public IAlgTool {
 public:
  static const InterfaceID& interfaceID() { return IID_INeighborXtalkTool; }

  /** \brief build internal map of Cal readout signals for use in subsequent calcXtalk() calls
      \note buildSignalMap() should be called @ beginning of each event.
  */
  virtual StatusCode buildSignalMap(const Event::CalDigiCol &digiCol) = 0;

  
  
  /** \brief calcuate total crosstalk for given cal diode from all neighboring crystals.

  Response is based on neigbhoring xtal signal provided in last call to BuildSignalMap()
  
  \param diodeIdx selected cal diode
  
  \return sum of xtalk for all neighboring xtals, 0 if not applicable.
                
  */
  virtual StatusCode calcXtalkCIDAC(CalUtil::DiodeIdx diodeIdx, float &xtalkDAC) const = 0;


  /** \brief calcuate total crosstalk for given cal diode from all neighboring crystals in units of MeV deposited @ center of xtal.

  Response is based on neigbhoring xtal signal provided in last call to BuildSignalMap()
  
  \param diodeIdx selected cal diode
  
  \return sum of xtalk for all neighboring xtals, 0 if not applicable.
                
  */
  virtual StatusCode calcXtalkMeV(CalUtil::DiodeIdx diodeIdx, float &xtalkMev) const = 0;


};

#endif
