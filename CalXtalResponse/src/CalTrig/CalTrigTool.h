//  $Header$
#ifndef CalTrigTool_h
#define CalTrigTool_h

// LOCAL
#include "CalXtalResponse/ICalTrigTool.h"
#include "CalXtalResponse/ICalCalibSvc.h"
#include "../CalCalib/IPrecalcCalibTool.h"
#include "CalXtalResponse/ICalSignalTool.h"

// GLAST
#include "CalUtil/CalDefs.h"
#include "CalUtil/CalArray.h"


// EXTLIB
#include "GaudiKernel/AlgTool.h"

// STD
#include <cstring>

class ICalSignalTool;
class IGlastDetSvc;

namespace Event {
  class GltDigi;
}


/*! \class CalTrigTool
  \author Zachary Fewtrell
  \brief Default implementation of ICalTrigTool.  Simulates GLAST Cal trigger
  response based on Cal xtal digi response or signal level

  provides multiple methods of calculating Cal or single crystal Trigger response.

  jobOptions:
  - CalCalibSvc       - cal calibration data source (default="CalCalibSvc")
  - PrecalcCalibTool  - source for derived calibration quantities (default="PrecalcCalibTool")
  - CalSignalToolName - source for cal signal levels (default="CalSignalTool")


*/

class CalTrigTool : public AlgTool, virtual public ICalTrigTool {
 public:
  /// default ctor, declares jobOptions
  CalTrigTool( const string& type, 
               const string& name, 
               const IInterface* parent);

  /// gets needed parameters and pointers to required services
  StatusCode initialize();

  StatusCode finalize() {return StatusCode::SUCCESS;}

  StatusCode calcXtalTrig(CalUtil::XtalIdx xtalIdx,
                          const CalUtil::CalArray<CalUtil::XtalRng, float> &adcPed,
                          CalUtil::CalArray<CalUtil::XtalDiode, bool> &trigBits,
                          Event::GltDigi *glt);

  StatusCode calcXtalTrig(const Event::CalDigi& calDigi,
                          CalUtil::CalArray<CalUtil::XtalDiode, bool> &trigBits,
                          Event::GltDigi *glt);


  StatusCode calcXtalTrig(CalUtil::XtalIdx xtalIdx,
                          const Event::CalDigi::CalXtalReadout &ro,
                          CalUtil::CalArray<CalUtil::XtalDiode, bool> &trigBits,
                          Event::GltDigi *glt);

  /// calculate trigger response for single crystal from CIDAC diode levels
  StatusCode calcXtalTrig(CalUtil::XtalIdx xtalIdx,
                          const ICalSignalTool::XtalSignalMap &cidac,
                          CalUtil::CalArray<CalUtil::XtalDiode, bool> &trigBits,
                          Event::GltDigi *glt
                          );


  StatusCode calcGlobalTrig(CalUtil::CalArray<CalUtil::DiodeNum, bool> &trigBits,
                            Event::GltDigi *glt);


  /// register GltDigi object in TDS if it has not been registered already.
  Event::GltDigi* setupGltDigi();

 private:

  /// calc full cal trigger response based on CalDigi
  StatusCode calcGlobalTrig(Event::CalDigiCol const &calDigiCol,
                            CalUtil::CalArray<CalUtil::DiodeNum, bool> &trigBits,
                            Event::GltDigi *glt);

  /// calc full cal trigger response based on McHits, using CalSignalTool
  StatusCode calcGlobalTrigSignalTool(CalUtil::CalArray<CalUtil::DiodeNum, bool> &trigBits,
                                      Event::GltDigi *glt);


  /// name of CalCalibSvc to use for calib constants.
  StringProperty  m_calCalibSvcName;      
  
  /// pointer to CalCalibSvc object.
  ICalCalibSvc   *m_calCalibSvc; 

  /// name of precalc calib tool
  StringProperty m_precalcCalibName;

  IPrecalcCalibTool *m_precalcCalibTool;

  /// ptr to event svc
  IDataProviderSvc* m_evtSvc;

  // name of CalSignalTool tool
  StringProperty m_calSignalToolName;
  /// ptr to CalSignalTool tool
  ICalSignalTool  *m_calSignalTool;

  /// used for constants & conversion routines.
  IGlastDetSvc* m_detSvc;

  /// list of active tower bays, populated at run time
  std::vector<CalUtil::TwrNum> m_twrList;

};


#endif
