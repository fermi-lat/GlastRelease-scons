#ifndef _TestAdcTool_H
#define _TestAdcTool_H 1

// Include files
#include "GaudiKernel/AlgTool.h"

#include "CalUtil/ICalAdcTool.h"
#include "CalUtil/ICalCalibSvc.h"
#include "CalUtil/IConvertAdc.h"

/*! \class TestAdcTool
  \author Zachary Fewtrell
  \brief Simple implementation of ICalAdcTool.  Faithfully pasted from CalDigiALg v1r3p7 and CalUtil/LinearConvertAdc v1r2p2
*/

class TestAdcTool : public AlgTool, virtual public ICalAdcTool {
public:
  /// default ctor, declares jobOptions
  TestAdcTool::TestAdcTool( const std::string& type, 
                            const std::string& name, 
                            const IInterface* parent);

  /// retrieves needed parameters and pointers to required services
  virtual StatusCode initialize();

  /// calculate xtal Adc response for all ranges basedon collection of McHits in xtal & diode regions
  StatusCode calculate(const idents::CalXtalId &xtalId, 
                       const std::vector<const Event::McIntegratingHit*> &hitList, 
                       bool &lacP,
                       bool &lacN,
                       idents::CalXtalId::AdcRange &rangeP, // output - best range
                       idents::CalXtalId::AdcRange &rangeN, // output - best range
                       std::vector<int> &adcP,              // output - ADC's for all ranges 0-3
                       std::vector<int> &adcN               // output - ADC's for all ranges 0-3
                       );
private:
  /// calculate light taper for the two xtal ends
  void calculateSignals(double relativePosition, 
                        double depositedEnergy,
                        std::pair<double, double> &signals);

  // *** STOLEN FROM CalDigi::CalDigiAlg *** //
  int m_xNum;                            ///< x tower number
  int m_eTowerCal;                       ///< detModel identifier for CAL
  int m_eLatTowers;                      ///< detModel identifier for LAT Towers
  double m_CsILength;                    ///< Xtal length
  int m_eDiodeMSmall;                    ///< detModel identifier for small minus-side diode
  int m_eDiodePSmall;                    ///< detModel identifier for small plus-side diode
  int m_eDiodeMLarge;                    ///< detModel identifier for large minus-side diode
  int m_eDiodePLarge;                    ///< detModel identifier for large plus-side diode
  int m_ePerMeV[2];                      ///< gain - electrons/MeV 1=Small, 0=Large
  double	m_ePerMeVinDiode;               ///< electrons per MeV for diodes
  int m_eXtal;
  double m_lightAtt;                     ///< light attenuation factor
  int m_maxAdc;                          ///< max value for ADC
  int m_nCsISeg;                         ///< number of geometric segments per Xtal
  int m_noise[2];                        ///< noise for diodes 1=Small, 0=Large units=electrons
  int m_pedestal;                        ///< single pedestal
  double m_thresh;                       ///< zero suppression threshold

  // *** STOLEN FROM CalUtil::LinearConvertAdc *** //
  double m_gain[2][4];                   ///< overall gain factor: MeV/channel
  double m_maxResponse[4];               ///< max energy responses per gain range

  /// names for volume identifier fields
  enum {fLATObjects, fTowerY, fTowerX, fTowerObjects, fLayer,
        fMeasure, fCALXtal,fCellCmp, fSegment};

  BooleanProperty m_doFluctuations;      ///< string flag for applying electron statistics fluctuations per channel
};

#endif //_TestAdcTool_H
