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

/*! @class TestEnergyTool
  \author Zachary Fewtrell
  \brief Simple implementation of ICalEnergyTool.  Faithfully pasted from CalXtalRecAlg v5r6p1
*/

class TestEnergyTool : public AlgTool, virtual public ICalEnergyTool {
public:

  /// default ctor, declares jobOptions.
  TestEnergyTool::TestEnergyTool( const std::string& type, 
                                  const std::string& name, 
                                  const IInterface* parent);


  /// retrieves needed parameters and pointers to required services
  virtual StatusCode initialize();

  /// calculate energy deposition given the digi response for both xtal faces
  StatusCode calculate(const idents::CalXtalId &xtalId, 
                       idents::CalXtalId::AdcRange rangeP,
                       idents::CalXtalId::AdcRange rangeN,
                       int adcP, 
                       int adcN,
                       float position,
                       float &energy                    // output
                       );
  
  /// calculate energy deposition given the digi response for one xtal face
  StatusCode calculate(const idents::CalXtalId &xtalId,
                       idents::CalXtalId::XtalFace face,
                       idents::CalXtalId::AdcRange range,
                       int adc, 
                       float position,
                       float &energy                    // output
                       );
private:
  int m_maxAdc;                          ///< max value for ADC
  double m_maxEnergy[4];                 ///< highest energy for each energy range
  int m_pedestal;                        ///< single pedestal
  int m_thresh;                          ///< zero suppression threshold

  CalibData::CalCalibPed* pPeds;         ///< pointer to pedestal data from CalCalibSvc
  CalibData::CalCalibGain* pGains;       ///< pointer to gain data from CalCalibSvc
  CalibData::CalCalibMuSlope* pMuSlopes; ///< pointer to muon slope data from CalCalibSvc

  /// Absolute time of first event (yyyy-mm-dd_hh:mm, trailing fields optional)
  StringProperty m_startTimeAsc;
  long m_startTime;                      ///< Absolute time of first event (seconds)
    
  std::string m_calibFlavor;             ///< "flavor" of calibration files
  IGlastDetSvc* detSvc;                  ///< pointer to the Glast Detector Service
  IDataProviderSvc* m_pCalibDataSvc;     ///< pointer to CalibDataSvc

  IDetDataSvc* m_detDataSvc;             ///< Handle to the IDetDataSvc interface of the CalibDataSvc

};

#endif //_TestEnergyTool_H
