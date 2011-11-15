#ifndef XtalDigiTool_h
#define XtalDigiTool_h
//  $Header$
// @file
//
//
// Author: Z.Fewtrell

// LOCAL
#include "CalXtalResponse/IXtalDigiTool.h"
#include "CalXtalResponse/ICalSignalTool.h"

// GLAST
#include "CalUtil/CalDefs.h"
#include "CalUtil/CalVec.h"
#include "CalUtil/CalArray.h"
#include "CalUtil/CalConfig.h"

// EXTLIB
#include "GaudiKernel/AlgTool.h"

// STD

// forward declarations
class ICalCalibSvc;
class TTree;
class ICalFailureModeSvc;
class IDataProviderSvc;
class IPrecalcCalibTool;

/*! \class XtalDigiTool
  \author Z.Fewtrell
  \brief Default implementation of IXtalDigiTool.  

  Convert Diode signals to CalDigi object for one Cal crystal.

  jobOptions:
  - CalCalibSvc (default="CalCalibSvc") - Cal calibration data source
  - PrecalcCalibTool (deefault="PrecalcCalibTool") - Derived Cal calibration data source.
  - CalSignalTool (default="CalSignalTool") - used to convert McIntegratingHits to diode signals
  
*/

class XtalDigiTool : public AlgTool, 
                     virtual public IXtalDigiTool
{
public:
  /// default ctor, declares jobOptions
  XtalDigiTool( const std::string& type, 
                const std::string& name, 
                const IInterface* parent);

  /// gets needed parameters and pointers to required services
  StatusCode initialize();

  StatusCode finalize();

  StatusCode calculate(Event::CalDigi &calDigi,
                       CalUtil::CalVec<CalUtil::FaceNum, bool> &lacBits,
                       bool zeroSuppress, string calFirstRng
                       );
private:
  
  /// select best adc range for both faces
  StatusCode XtalDigiTool::rangeSelect(const CalUtil::XtalIdx xtalIdx,
                                       CalUtil::CalArray<CalUtil::XtalRng, float> &adcPed,
                                       CalUtil::CalVec<CalUtil::FaceNum, CalUtil::RngNum> &bestRng);

  /// get failureMode bits
  StatusCode getFailureStatus(const idents::CalXtalId xtalId,
                              const CalUtil::CalVec<CalUtil::FaceNum, bool> &lacBits,
                              unsigned short &failureStatus
                              );

  /// populate Digi TDS class
  StatusCode fillDigi(Event::CalDigi &calDigi,
                      const CalUtil::CalArray<CalUtil::XtalRng, float> &adcPed,
                      const CalUtil::CalVec<CalUtil::FaceNum, CalUtil::RngNum> &bestRng,
                      unsigned short failureStatus
                      );

  /// name of CalCalibSvc to use for calib constants.
  StringProperty  m_calCalibSvcName;      
  
  /// pointer to CalCalibSvc object.
  ICalCalibSvc   *m_calCalibSvc; 

  /// max val for ADC
  int m_maxAdc;  

  /// pointer to optional CalFailureMode
  ICalFailureModeSvc *m_calFailureModeSvc;

  /// name of precalc calib tool
  StringProperty m_precalcCalibName;
  
  /// pointer to precalcCalibTool
  IPrecalcCalibTool *m_precalcCalib;

  /// name of CalSignalTool tool, used for calculating diode signal levels from McIntegratingHits
  StringProperty m_calSignalToolName;

  /// ptr to CalSignalTool tool
  ICalSignalTool  *m_calSignalTool;


};


#endif
