#ifndef ICalDiagnosticTool_H
#define ICalDiagnosticTool_H
// $Header$
/** @file
    @author Z.Fewtrell
*/
/*! @class ICalDiagnosticTool
 *
 * \author Z.Fewtrell
 *
 * \brief Abstract interface class for determination of Cal Diagnostic
 * trigger primitives
 *
 *
 */

// GLAST INCLUDES
#include "CalUtil/CalDefs.h"
#include "LdfEvent/DiagnosticData.h"

// EXTLIB INCLUDES
#include "GaudiKernel/IAlgTool.h"

// STD INCLUDES
#include <memory>

static const InterfaceID IID_ICalDiagnosticTool("ICalDiagnosticTool", 1, 0);

class ICalDiagnosticTool : virtual public IAlgTool {
public:

  static const InterfaceID& interfaceID() { return IID_ICalDiagnosticTool; }

  virtual ~ICalDiagnosticTool() {};

  /// retrieve CalDiagnosticData object for given Tower bay and Cal
  /// layer
  /// @return Null pointer on error.
  virtual std::auto_ptr<LdfEvent::CalDiagnosticData> getDiagnosticData(const CalUtil::TwrNum twr,
                                                                       const CalUtil::LyrNum lyr) =0;

};

#endif //
