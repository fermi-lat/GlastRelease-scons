#ifndef _TestPosTool_H
#define _TestPosTool_H 1

// Include files
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/IDetDataSvc.h"

#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "CalibData/Cal/CalCalibMuSlope.h"

#include "CalUtil/ICalPosTool.h"
#include "CalUtil/ICalEnergyTool.h"
#include "CalUtil/ICalCalibSvc.h"

class TestPosTool : public AlgTool, virtual public ICalPosTool {
public:
  TestPosTool::TestPosTool( const std::string& type, 
                            const std::string& name, 
                            const IInterface* parent);

  virtual StatusCode initialize();


  // calculate position given the digital response on both faces
  StatusCode calculate(const idents::CalXtalId &xtalId,
                       int adcP, 
                       idents::CalXtalId::AdcRange rangeP,
                       int adcN, 
                       idents::CalXtalId::AdcRange rangeN,
                       float &position                // output
                       );
private:
  IGlastDetSvc* detSvc; ///< pointer to the Glast Detector Service
  /// constants defining the position of the fields in VolumeIdentifier 
  enum {fLATObjects, fTowerY, fTowerX, fTowerObjects, fLayer,
        fMeasure, fCALXtal,fCellCmp, fSegment};

  double m_CsILength;  ///< Xtal length
  double m_CsIWidth;  ///< Xtal width
  double m_CsIHeight;  ///< Xtal height
  int m_eTowerCAL;
  ///< the value of fTowerObject field, defining calorimeter module 
  int m_eLATTowers; ///< the value of fLATObjects field, defining LAT towers 
  int m_eXtal;      ///< the value of fCellCmp field defining CsI crystal
  double m_lightAtt;  ///< light attenuation factor
  int m_nCsISeg;  ///< number of geometric segments per Xtal
  int m_xNum;  ///< x tower number
  int m_yNum;  ///< y tower number

  CalibData::CalCalibMuSlope* pMuSlopes;

  /// "flavor" of calibration files
  std::string m_calibFlavor;

  ICalEnergyTool *m_pCalEnergyTool;

  /// Absolute time of first event (yyyy-mm-dd_hh:mm, trailing fields
  /// optional)
  std::string m_startTimeAsc;

  /// Absolute time of first event (seconds)
  long m_startTime;

  IDataProviderSvc* m_pCalibDataSvc;

  /// Handle to the IDetDataSvc interface of the CalibDataSvc
  IDetDataSvc* m_detDataSvc;


};

#endif //_TestPosTool_H
