//  $Header$
/** @file
    @author Z.Fewtrell
*/

#ifndef CalDiagnosticTool_h
#define CalDiagnosticTool_h

// LOCAL
#include "CalXtalResponse/ICalDiagnosticTool.h"

// GLAST
#include "CalUtil/CalDefs.h"
#include "CalUtil/CalVec.h"
#include "CalUtil/CalDiagnosticWord.h"

// EXTLIB
#include "GaudiKernel/AlgTool.h"


// STD

class IPrecalcCalibTool;
class ICalSignalTool;
class ICalTrigTool;

/*! \class CalDiagnosticTool
  \author Z.Fewtrell
  \brief Default implementation of ICalDiagnosticTool.  Generate Cal
  Diagnostic trigger primitives from CalTrigTool, CalSignalTool output
  & current calibrations.

  jobOptions:
  - PrecalcCalibTool  - source for derived calibration quantities, used for LAC thresholds (default="PrecalcCalibTool")
  - CalTrigToolName   - source for cal trigger response (default="CalTrigTool"
  - CalSignalToolName - source for cal signal levels, used for LAC bits  (default="CalSignalTool")


*/

class CalDiagnosticTool : 
  public AlgTool, 
  virtual public ICalDiagnosticTool
{
public:
  /// default ctor, declares jobOptions
  CalDiagnosticTool( const std::string& type, 
                     const std::string& name, 
                     const IInterface* parent);

  /// gets needed parameters and pointers to required services
  StatusCode initialize();

  StatusCode finalize() {return StatusCode::SUCCESS;}

  /// retrieve CalDiagnosticData object for given Tower bay and Cal
  /// layer
  /// @return Null pointer on error.
  ///
  /// create  new CalDiagnosticData object for each call, fill bits
  /// from information in CalTrigTool, CalSignalTool & PrecalcCalib
  virtual std::auto_ptr<LdfEvent::CalDiagnosticData> getDiagnosticData(const CalUtil::TwrNum twr,
                                                                       const CalUtil::LyrNum lyr);
private:
  /// calculate LAC bits for single layer
  StatusCode calcLACBits(const CalUtil::TwrNum twr,
                         const CalUtil::LyrNum lyr,
                         CalUtil::CalDiagnosticWord::CalDiagLACBits &lacBits);

  /// calculate OR'd trigger bits for single cal diagnostic 
  StatusCode calcTrigBits(const CalUtil::TwrNum twr,
                          const CalUtil::LyrNum lyr,
                          CalUtil::CalDiagnosticWord::CalDiagTrigBits &trigBits);

  /// name of precalc calib tool
  StringProperty m_precalcCalibName;

  IPrecalcCalibTool *m_precalcCalibTool;

  // name of CalSignalTool tool
  StringProperty m_calSignalToolName;
  /// ptr to CalSignalTool tool
  ICalSignalTool  *m_calSignalTool;

  StringProperty m_calTrigToolName;
  ICalTrigTool *m_calTrigTool;

};


#endif
