#ifndef _TestAdcTool_H
#define _TestAdcTool_H 1

// Include files
#include "GaudiKernel/AlgTool.h"

#include "CalUtil/ICalAdcTool.h"
#include "CalUtil/ICalCalibSvc.h"
#include "CalUtil/IConvertAdc.h"

class TestAdcTool : public AlgTool, virtual public ICalAdcTool {
public:
  TestAdcTool::TestAdcTool( const std::string& type, 
                            const std::string& name, 
                            const IInterface* parent);

  virtual StatusCode initialize();

  StatusCode calculate(const idents::CalXtalId &xtalId, 
                       const std::vector<const Event::McIntegratingHit*> &hitList,
                       std::vector<int> &adcP,              // output - ADC's for all ranges 0-3
                       idents::CalXtalId::AdcRange &rangeP, // output - best range
                       std::vector<int> &adcN,              // output - ADC's for all ranges 0-3
                       idents::CalXtalId::AdcRange &rangeN  // output - best range
                       );
private:
   
  void calculateSignals(double relativePosition, 
                        double depositedEnergy,
                        std::pair<double, double> &signals);

  // *** STOLEN FROM CalDigi::CalDigiAlg *** //
  /// Xtal length
  double m_CsILength; 
  /// detModel identifier for small minus-side diode
  int m_eDiodeMSmall;   
  /// detModel identifier for small plus-side diode
  int m_eDiodePSmall;
  /// detModel identifier for large minus-side diode
  int m_eDiodeMLarge;
  /// detModel identifier for large plus-side diode
  int m_eDiodePLarge;
  /// gain - electrons/MeV 1=Small, 0=Large
  int m_ePerMeV[2];  
  /// electrons per MeV for diodes
  double	m_ePerMeVinDiode;
  int m_eXtal;
  /// light attenuation factor
  double m_lightAtt;  
  /// max value for ADC
  int m_maxAdc;  
  /// number of geometric segments per Xtal
  int m_nCsISeg;  
  /// noise for diodes 1=Small, 0=Large units=electrons
  int m_noise[2];  
  /// single pedestal
  int m_pedestal;  

  // *** STOLEN FROM CalUtil::LinearConvertAdc *** //
  /// overall gain factor: MeV/channel
  double m_gain[2][4];
  /// max energy responses per gain range
  double m_maxResponse[4];

  /// names for identifier fields
  enum {fLATObjects, fTowerY, fTowerX, fTowerObjects, fLayer,
        fMeasure, fCALXtal,fCellCmp, fSegment};

  /// string flag for applying electron statistics fluctuations per channel
  BooleanProperty m_doFluctuations;
  /// input XML file containing parameters for Digitization
  StringProperty m_xmlFile;

};

#endif //_TestAdcTool_H
