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

/*! @class TestPosTool
  \author Zachary Fewtrell
  \brief Simple implementation of ICalPosTool.  Faithfully pasted from CalXtalRecAlg v5r6p1
*/
class TestPosTool : public AlgTool, virtual public ICalPosTool {
public:
  
  /// default ctor, declares jobOptions
  TestPosTool::TestPosTool( const std::string& type, 
                            const std::string& name, 
                            const IInterface* parent);

  /// retrieves needed paramters and pointers to required services
  virtual StatusCode initialize();


  /// calculate position of energy deposition relative to xtal-center = 0
  StatusCode calculate(const idents::CalXtalId &xtalId,
                       int adcP, 
                       idents::CalXtalId::AdcRange rangeP,
                       int adcN, 
                       idents::CalXtalId::AdcRange rangeN,
                       float &position                // output
                       );
private:
  IGlastDetSvc* detSvc;                  ///< pointer to the Glast Detector Service

  /// constants defining the position of the fields in VolumeIdentifier 
  enum {fLATObjects, fTowerY, fTowerX, fTowerObjects, fLayer,
        fMeasure, fCALXtal,fCellCmp, fSegment};

  double m_CsILength;                    ///< Xtal length
  double m_CsIWidth;                     ///< Xtal width
  double m_CsIHeight;                    ///< Xtal height
  int m_eTowerCAL;
                                         ///< the value of fTowerObject field, defining calorimeter module 
  int m_eLATTowers;                      ///< the value of fLATObjects field, defining LAT towers 
  int m_eXtal;                           ///< the value of fCellCmp field defining CsI crystal
  double m_lightAtt;                     ///< light attenuation factor
  int m_nCsISeg;                         ///< number of geometric segments per Xtal
  int m_xNum;                            ///< x tower number
  int m_yNum;                            ///< y tower number

  CalibData::CalCalibMuSlope* pMuSlopes; ///< pointer to muon slope calibration data in TDS

  std::string m_calibFlavor;             ///< "flavor" of calibration files


  ICalEnergyTool *m_pCalEnergyTool;      ///< CalEnergyTool is needed to obtain energy estimate for calculation

  /// Absolute time of first event (yyyy-mm-dd_hh:mm, trailing fields optional)
  StringProperty m_startTimeAsc;
  long m_startTime;                      ///< Absolute time of first event (seconds)


  IDataProviderSvc* m_pCalibDataSvc;     ///< pointer to CalibDataSvc, where i get my calib constants
  IDetDataSvc* m_detDataSvc;             ///< Handle to the IDetDataSvc interface of the CalibDataSvc
};

#endif //_TestPosTool_H
