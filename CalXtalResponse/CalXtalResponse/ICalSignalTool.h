#ifndef ICalSignalTool_H
#define ICalSignalTool_H
//  $Header$

// LOCAL INCLUDES


// GLAST INCLUDES
#include "CalUtil/CalDefs.h"
#include "CalUtil/CalArray.h"
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
static const InterfaceID IID_ICalSignalTool("ICalSignalTool", 0, 1);

/*! @class ICalSignalTool
 * \brief Abstract interface for calculation of Cal diode signal levels, also maitain relational map between crystals and MCIntegratingHits
 * \author Zach Fewtrell
 *
 */
class ICalSignalTool : virtual public IInterface {
 public:
  static const InterfaceID& interfaceID() { return IID_ICalSignalTool; }

  /// map used to associate Cal CsI crystals with electronic signal
  typedef CalUtil::CalVec<CalUtil::DiodeIdx, float> CalSignalMap;


  /// \brief return map of diodes and associated signal levels
  /// \note map is only valid until end of current event
  /// \return null on error
  virtual const CalSignalMap *getCalSignalMap() = 0;

  /// contains all signals for single crystal
  typedef CalUtil::CalArray<CalUtil::XtalDiode, float> XtalSignalMap;

  /// \brief return all diode signals for single xtal
  /// \param xtalIdx selected Cal crytsal
  /// \param output signal levels for each xtal diode
  virtual StatusCode getXtalSignalMap(CalUtil::XtalIdx xtalIdx,
                                      XtalSignalMap &xtalSignalMap) = 0;

  /// map used to associate Cal MCIntegrating hits with crystals
  typedef std::multimap<CalUtil::XtalIdx, Event::McIntegratingHit* > CalRelationMap;
  typedef CalRelationMap::value_type CalRelation;

  /// \brief return map of Cal Crystals and associated McIntegrating Hits
  /// \notes crystals with 0 hits are omitted.
  /// \note map is only valid until end of current event
  /// \return null on error
  virtual const CalRelationMap *getCalRelationMap() = 0;

protected:
  friend class test_CalXtalResponse;

  /// clear internal tables (used by test app, unneeded during normal operations)
  virtual void newEvent() = 0;

};

#endif // ICalSignalTool_H
