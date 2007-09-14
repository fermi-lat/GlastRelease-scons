//  $Header$
#ifndef CalTrigTool_h
#define CalTrigTool_h

// LOCAL
#include "CalXtalResponse/ICalTrigTool.h"
#include "CalXtalResponse/ICalCalibSvc.h"
#include "../CalCalib/IPrecalcCalibTool.h"

// GLAST
#include "CalUtil/CalDefs.h"
#include "CalUtil/CalArray.h"


// EXTLIB
#include "GaudiKernel/AlgTool.h"

// STD
#include <cstring>

/*! \class CalTrigTool
  \author Zachary Fewtrell
  \brief Official implementation of ICalTrigTool.  Simulates GLAST Cal trigger
  response based on Cal xtal digi response.


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

  StatusCode calcXtalTrig(CalUtil::XtalIdx xtalIdx,
                          const CalUtil::CalArray<CalUtil::XtalDiode, float> &cidac,
                          CalUtil::CalArray<CalUtil::XtalDiode, bool> &trigBits,
                          Event::GltDigi *glt
                          );


  StatusCode calcGlobalTrig(const Event::CalDigiCol& calDigiCol,
                            CalUtil::CalArray<CalUtil::DiodeNum, bool> &trigBits,
                            Event::GltDigi *glt);


  Event::GltDigi* setupGltDigi(IDataProviderSvc *eventSvc);

 private:

  /// name of CalCalibSvc to use for calib constants.
  StringProperty  m_calCalibSvcName;      
  
  /// pointer to CalCalibSvc object.
  ICalCalibSvc   *m_calCalibSvc; 

  /// name of precalc calib tool
  StringProperty m_precalcCalibName;

  IPrecalcCalibTool *m_precalcCalibTool;

};


#endif
