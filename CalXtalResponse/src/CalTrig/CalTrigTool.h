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

  StatusCode calcXtalTrig(XtalIdx xtalIdx,
                          const CalArray<XtalRng, float> &adcPed,
                          CalArray<XtalDiode, bool> &trigBits,
                          Event::GltDigi *glt);

  StatusCode calcXtalTrig(const Event::CalDigi& calDigi,
                          CalArray<XtalDiode, bool> &trigBits,
                          Event::GltDigi *glt);


  StatusCode calcXtalTrig(XtalIdx xtalIdx,
                          const Event::CalDigi::CalXtalReadout &ro,
                          CalArray<XtalDiode, bool> &trigBits,
                          Event::GltDigi *glt);

  StatusCode calcXtalTrig(XtalIdx xtalIdx,
                          const CalArray<XtalDiode, float> &cidac,
                          CalArray<XtalDiode, bool> &trigBits,
                          Event::GltDigi *glt
                          );


  StatusCode calcGlobalTrig(const Event::CalDigiCol& calDigiCol,
                            CalArray<DiodeNum,bool> &trigBits,
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
