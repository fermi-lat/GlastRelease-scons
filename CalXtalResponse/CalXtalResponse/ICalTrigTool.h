#ifndef ICalTrigTool_H
#define ICalTrigTool_H
// $Header$
/** @file
    @author Z.Fewtrell
*/
/*! @class ICalTrigTool
 *
 * \author Z.Fewtrell
 *
 * \brief Abstract interface class for determination of Cal FLE & FHE trigger response 
 *
 *
 */

// GLAST INCLUDES
#include "Event/TopLevel/Event.h"
#include "CalUtil/CalDefs.h"
#include "CalUtil/CalVec.h"

// EXTLIB INCLUDES
#include "GaudiKernel/IAlgTool.h"

static const InterfaceID IID_ICalTrigTool("ICalTrigTool", 2, 0);

class ICalTrigTool : virtual public IAlgTool {
public:

  static const InterfaceID& interfaceID() { return IID_ICalTrigTool; }

  virtual ~ICalTrigTool() {};

  /// return trigger response for given channel (specify xtal, face & diode)
  virtual StatusCode getTriggerBit(CalUtil::DiodeIdx diodeIdx, bool &trigBit) = 0;

  /// \brief return 16 bit trigger vector for given diode size, one bit per tower
  virtual StatusCode getCALTriggerVector(idents::CalXtalId::DiodeType diode, 
                                         unsigned short &vec) = 0;

};

#endif //
