#ifndef _TestEnergyTool_H
#define _TestEnergyTool_H 1

// Include files
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/IDataProviderSvc.h"

#include "CalibData/Cal/CalCalibPed.h"
#include "CalibData/Cal/CalCalibGain.h"
#include "CalibData/Cal/CalCalibMuSlope.h"
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"

#include "CalUtil/ICalEnergyTool.h"
#include "CalUtil/ICalCalibSvc.h"

class TestEnergyTool : public AlgTool, virtual public ICalEnergyTool {
public:
  TestEnergyTool::TestEnergyTool( const std::string& type, 
                                  const std::string& name, 
                                  const IInterface* parent);


  virtual StatusCode initialize();

  StatusCode calculate(const idents::CalXtalId &xtalId, 
                       idents::CalXtalId::AdcRange rangeP,
                       idents::CalXtalId::AdcRange rangeN,
                       int adcP, 
                       int adcN,
                       float position,
                       float &energy                    // output
                       );
  
  // calculate energy from xtalId, one face/range/adc, and a position
  StatusCode calculate(const idents::CalXtalId &xtalId,
                       idents::CalXtalId::AdcRange range,
                       idents::CalXtalId::XtalFace face,
                       int adc, 
                       float position,
                       float &energy                    // output
                       );
private:
  int m_maxAdc;  ///< max value for ADC
  double m_maxEnergy[4];  ///< highest energy for each energy range
  int m_pedestal;  ///< single pedestal
  int m_thresh;  ///< zero suppression threshold

  CalibData::CalCalibPed* pPeds;
  CalibData::CalCalibGain* pGains;
  CalibData::CalCalibMuSlope* pMuSlopes;

  /// Absolute time of first event (yyyy-mm-dd_hh:mm, trailing fields
  /// optional)
  std::string m_startTimeAsc;

  /// Absolute time of first event (seconds)
  long m_startTime;
    
  /// "flavor" of calibration files
  std::string m_calibFlavor;

  IGlastDetSvc* detSvc; ///< pointer to the Glast Detector Service

  IDataProviderSvc* m_pCalibDataSvc;

  /// Handle to the IDetDataSvc interface of the CalibDataSvc
  IDetDataSvc* m_detDataSvc;
};

#endif //_TestEnergyTool_H
