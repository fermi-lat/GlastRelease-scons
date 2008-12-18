#ifndef ICalSignalTool_H
#define ICalSignalTool_H
//  $Header$
/** @file
    @author Z.Fewtrell
*/

// LOCAL INCLUDES


// GLAST INCLUDES
#include "CalUtil/CalDefs.h"
#include "CalUtil/CalVec.h"

// EXTLIB INCLUDES
#include "GaudiKernel/IInterface.h"

// STD INCLUDES
#include <map>

namespace Event {
  class McIntegratingHit;
}

class test_CalXtalResponse;


// Declaration of the interface ID ( interface id, major version,
// minor version)
static const InterfaceID IID_ICalSignalTool("ICalSignalTool", 0, 2);

/*! @class ICalSignalTool
 * \brief Abstract interface for calculation of Cal diode signal
 * levels from MC Hits, also maintain relational map between crystals and 
 * MCIntegratingHits
 * 
 * \author Z.Fewtrell
 *
 */
class ICalSignalTool : virtual public IInterface {
 public:
  static const InterfaceID& interfaceID() { return IID_ICalSignalTool; }

  virtual ~ICalSignalTool() {};

  /// \brief return signal in charge-injection DAC units for single
  /// crystal diode.
  virtual StatusCode getDiodeSignal(CalUtil::DiodeIdx, float &signal) = 0;

  /// \brief return signal in charge-injection DAC units for single
  /// crystal diode.
  virtual StatusCode getTrigDiodeSignal(CalUtil::DiodeIdx, float &signal) = 0;

  /// map used to associate Cal MCIntegrating hits with crystals
  typedef std::multimap<CalUtil::XtalIdx, Event::McIntegratingHit* > CalRelationMap;
  typedef CalRelationMap::value_type CalRelation;

  /// \brief return map of Cal Crystals and associated McIntegrating Hits
  /// \notes crystals with 0 hits are omitted.
  /// \note map is only valid until end of current event
  /// \return null on error
  virtual const CalRelationMap *getCalRelationMap() = 0;
};

#endif // ICalSignalTool_H
