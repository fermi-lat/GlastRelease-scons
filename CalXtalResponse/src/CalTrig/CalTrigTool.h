//  $Header$
/** @file
    @author Z.Fewtrell
*/

#ifndef CalTrigTool_h
#define CalTrigTool_h



// LOCAL
#include "CalXtalResponse/ICalTrigTool.h"
#include "CalXtalResponse/ICalCalibSvc.h"
#include "CalXtalResponse/ICalSignalTool.h"


// GLAST
#include "CalUtil/CalDefs.h"
#include "CalUtil/CalVec.h"
#include "Event/Digi/CalDigi.h"

// EXTLIB
#include "GaudiKernel/IIncidentListener.h"
#include "GaudiKernel/AlgTool.h"


// STD

class IGlastDetSvc;
class IPrecalcCalibTool;
class IDataProviderSvc;

/*! \class CalTrigTool
  \author Z.Fewtrell
  \brief Default implementation of ICalTrigTool.  Simulates GLAST Cal trigger
  response based on either Cal xtal digi response or simulated signal
  level (from MC)


  jobOptions:
  - CalCalibSvc       - cal calibration data source (default="CalCalibSvc")
  - PrecalcCalibTool  - source for derived calibration quantities (default="PrecalcCalibTool")
  - CalSignalToolName - source for cal signal levels (default="CalSignalTool")


*/

class CalTrigTool : 
  public AlgTool, 
  virtual public ICalTrigTool,
  virtual public IIncidentListener
{
public:
  /// default ctor, declares jobOptions
  CalTrigTool( const std::string& type, 
               const std::string& name, 
               const IInterface* parent);

  /// gets needed parameters and pointers to required services
  StatusCode initialize();

  StatusCode finalize() {return StatusCode::SUCCESS;}

  /// \brief return 16 bit trigger vector for FLE trigger, one bit per tower
  StatusCode getCALTriggerVector(idents::CalXtalId::DiodeType diode, 
                                 unsigned short &vec);

  /// return trigger response for given channel (specify xtal, face & diode)
  StatusCode getTriggerBit(CalUtil::DiodeIdx diodeIdx, bool &trigBit);

  /// hook the BeginEvent so that we can check our validity once per event.
  void handle ( const Incident& inc );

private:

  /// update all Cal Trigger bits w/ current TDS data (either MC or
  /// Digi)
  StatusCode calcGlobalTrig();

  /// update all Cal Trigger bits w/ current TDS MC data (via CalSignalTool)
  StatusCode CalTrigTool::calcGlobalTrigSignalTool();

  /// update all Cal Trigger bits w/ current TDS data from CalDigi
  StatusCode CalTrigTool::calcGlobalTrigDigi(const Event::CalDigiCol &calDigiCol);

  /// call to clear all trigger bits (like @ beginning of event)
  void newEvent();

  /// set single Cal trigger channel to high, set all corresponding
  /// records simultaneously.
  void setSingleBit(const CalUtil::DiodeIdx diodeIdx);

  /// calculate single xtal trigger response from cal Digi object.
  /// 
  /// store results in internal private tables
  StatusCode calcXtalTrig(const Event::CalDigi& calDigi);

  /// calculate trigger response for single crystal from CIDAC diode
  /// levels
  /// 
  /// store results in internal private tables
  StatusCode calcXtalTrigSignalTool(const CalUtil::XtalIdx xtalIdx);

  /// calculate trigger response from single crystal readout
  StatusCode calcXtalTrig(const CalUtil::XtalIdx xtalIdx,
                          const Event::CalDigi::CalXtalReadout &ro);

  /// calculate trigger response from pedestal subtracted adc values
  /// for all 8 adc channels
  StatusCode calcXtalTrig(const CalUtil::XtalIdx xtalIdx,
                          const CalUtil::CalVec<CalUtil::XtalRng, float> &adcPed);

  /// store trigger values for each channel in Cal
  typedef CalUtil::CalVec<CalUtil::DiodeIdx, bool> CalTriggerMap;

  /// store trigger values for each channel in Cal
  CalTriggerMap m_calTriggerMap;

  /// store FLE trigger vector for each tower in cal
  CalUtil::CalVec<CalUtil::DiodeNum,
                  unsigned short> m_calTriggerVec;

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

  /// if current private data store is valid
  bool m_isValid;
  
  /// CAL rows/columns selection for FLE/FHE trigger generation:
  /// "default" - enable all channels to produce trigger
  /// "erec"    - enable only even columns in even rows and odd columns in odd rows
  /// "eroc"    - enable only odd columns in even rows and even columns in odd rows
  StringProperty m_selectionRule;
};


#endif
