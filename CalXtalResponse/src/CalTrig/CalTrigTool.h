#ifndef CalTrigTool_h
#define CalTrigTool_h

// LOCAL
#include "CalXtalResponse/ICalTrigTool.h"
#include "CalXtalResponse/ICalCalibSvc.h"

// GLAST
#include "CalUtil/CalDefs.h"
#include "CalUtil/ICalFailureModeSvc.h"
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
  virtual StatusCode initialize();

  virtual StatusCode finalize();

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



  StatusCode calcGlobalTrig(const Event::CalDigiCol& calDigiCol,
                            CalArray<DiodeNum,bool> &trigBits,
                            Event::GltDigi *glt);


  Event::GltDigi* setupGltDigi(IDataProviderSvc *eventSvc);

 private:
  /// \brief Convert LEX8 adc value to LEX1 scale for given xtal face
  /// \note support LEX8 values > 4095 via extrapolation
  /// \param uldTholdX8 is needed for the calculation.  i know you already
  /// have it, so there's no sense in me retrieving it again.
  StatusCode lex8_to_lex1(FaceIdx faceIdx, 
                          float l8adc, float &l1adc);

  /// \brief Convert HEX8 adc value to HEX1 scale for given xtal face
  /// \note support HEX8 values > 4095 via extrapolation
  /// \param uldTholdX8 is needed for the calculation.  i know you already
  /// have it, so there's no sense in me retrieving it again.
  StatusCode hex8_to_hex1(FaceIdx faceIdx, 
                          float h8adc, float &h1adc);

  /// name of CalCalibSvc to use for calib constants.
  StringProperty  m_calCalibSvcName;      
  
  /// pointer to CalCalibSvc object.
  ICalCalibSvc   *m_calCalibSvc; 

  /// used to hold interim variables for current xtal calc
  struct AlgData {
    void Clear() {memset(this, 0, sizeof(AlgData));}

    TwrNum twr;
    LyrNum lyr;
    ColNum col;

    CalArray<XtalDiode, float> trigThresh;
    CalArray<XtalRng, float> uldThold;
    
    
  };

  /// used to hold interim variables for current xtal calc
  AlgData m_dat;

};


#endif
