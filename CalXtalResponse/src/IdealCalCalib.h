#ifndef IdealCalCalib_H
#define IdealCalCalib_H

// LOCAL INCLUDES

// GLAST INCLUDES

// EXTLIB
#include "GaudiKernel/StatusCode.h"

// STD
#include <vector>
#include <string>

using namespace std;

/**
   \brief Contains all values necessary to generate full set of ideal cal
   calibration constants

   contains code to read constants in from IFile XML file.

*/
class IdealCalCalib {
 public:
  /// read in calibration values from IFile XML file
  StatusCode readCfgFile(const string &path);

  //////////////////
  //// PEDS ///
  //////////////////
  vector<double> pedVals;       ///< 1 pedestal value per range
  vector<double> pedCos;        ///< associated correlated pedestal cosines
  double         pedSigPct;     ///< sigma/val for all pedestal calibrations

  //////////////////
  //// ASYMMETRY ///
  //////////////////
  double asymLrgNeg;            ///< val for Large Diode asym spline at neg face
  double asymLrgPos;            ///< val for Large Diode asym spline at pos face
        
  double asymSmNeg;             ///< val for Small Diode asym spline at neg face
  double asymSmPos;             ///< val for Small Diode asym spline at pos face
        
  double asymPSNBNeg;           ///< val for PSNB asym spline at neg face
  double asymPSNBPos;           ///< val for PSNB asym spline at pos face
        
  double asymNSPBNeg;           ///< val for NSPB asym spline at neg face
  double asymNSPBPos;           ///< val for NSPB asym spline at pos face
        
  double asymSigPct;            ///< sigma/val for all vals in asym group 
  

  ///////////////////
  //// MEVPERDAC ////
  ///////////////////
  double mpdLrg;                ///< Large diode MevPerDac 
  double mpdSm;                 ///< Small diode MevPerDac 
        
  double mpdSigPct;             ///< sigma/val for all vals in MPD group 
  

  /////////////////
  //// THOLD_CI ///
  /////////////////
  double ciFLE;                 ///< FLE threshold 
  double ciFHE;                 ///< FHE threshold 
  double ciLAC;                 ///< LAC threshold 
  vector<double> ciULD;         ///< 1 ULD threshold per ADC range
  double ciSigPct;              ///< sigma/val for all vals in THOLD_CI group
  vector<double> ciPeds;        ///< 1 pedestal value per ADC range           
  

  ///////////////////
  //// THOLD_MUON ///
  ///////////////////
  double muonFLE;               ///< FLE threshold 
  double muonFHE;               ///< FHE threshold 
  double muonSigPct;            ///< sigma/val for all vals in THOLD_MUON group 
  vector<double> muonPeds;      ///< 1 pedestal value per ADC range      

  ///////////////////
  //// INT_NONLIN ///
  ///////////////////
  vector<double> inlADCPerDAC;          ///< ADC/DAC for each ADC range 
  double inlSigPct;             ///< sigma/val for all vals in INT_NONLIN group

 private:

  //-- SECTION DECRIPTION STRINGS --//
  static const string PEDS;
  static const string ASYM;
  static const string THOLD_CI;
  static const string THOLD_MUON;
  static const string INL;
  static const string MPD;
};



#endif
