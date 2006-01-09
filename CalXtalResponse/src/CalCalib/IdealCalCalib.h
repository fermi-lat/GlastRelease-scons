#ifndef IdealCalCalib_H
#define IdealCalCalib_H
// $Header$

// LOCAL INCLUDES

// GLAST INCLUDES

// EXTLIB
#include "GaudiKernel/StatusCode.h"

// STD
#include <vector>
#include <string>

using namespace std;

/** @class IdealCalCalib
    @author Zachary Fewtrell
    \brief Contains all values necessary to generate full set of ideal cal
    calibration constants
    
    contains code to read constants in from IFile XML file.
*/
class IdealCalCalib {
 public:
  IdealCalCalib();

  /// read in calibration values from IFile XML file
  StatusCode readCfgFile(const string &path);

  //////////////////
  //// PEDS ///
  //////////////////
  /// 1 pedestal value per range
  vector<double> pedVals;       
  /// associated correlated pedestal cosines
  vector<double> pedCos;        
  /// sigma/val for all pedestal calibrations
  vector<double> pedSigs;       

  //////////////////
  //// ASYMMETRY ///
  //////////////////
  /// val for Large Diode asym spline at neg face
  float asymLrgNeg;            
  /// val for Large Diode asym spline at pos face
  float asymLrgPos;            
        
  /// val for Small Diode asym spline at neg face
  float asymSmNeg;             
  /// val for Small Diode asym spline at pos face
  float asymSmPos;             
        
  /// val for PSNB asym spline at neg face
  float asymPSNBNeg;           
  /// val for PSNB asym spline at pos face
  float asymPSNBPos;           
        
  /// val for NSPB asym spline at neg face
  float asymNSPBNeg;           
  /// val for NSPB asym spline at pos face
  float asymNSPBPos;           
        
  /// sigma/val for all vals in asym group 
  float asymSigPct;            
  

  ///////////////////
  //// MEVPERDAC ////
  ///////////////////
  /// Large diode MevPerDac 
  float mpdLrg;                
  /// Small diode MevPerDac 
  float mpdSm;                 
        
  /// sigma/val for all vals in MPD group 
  float mpdSigPct;             
  

  /////////////////
  //// THOLD_CI ///
  /////////////////
  /// FLE threshold 
  float ciFLE;                 
  /// FHE threshold 
  float ciFHE;                 
  /// LAC threshold 
  float ciLAC;                 
  /// 1 ULD threshold per ADC range
  vector<double> ciULD;         
  /// sigma/val for all vals in THOLD_CI group
  float ciSigPct;              
  /// 1 pedestal value per ADC range           
  vector<double> ciPeds;        
  

  ///////////////////
  //// THOLD_MUON ///
  ///////////////////
  /// FLE threshold 
  float muonFLE;               
  /// FHE threshold 
  float muonFHE;               
  /// sigma/val for all vals in THOLD_MUON group 
  float muonSigPct;            
  /// 1 pedestal value per ADC range      
  vector<double> muonPeds;      

  ///////////////////
  //// INT_NONLIN ///
  ///////////////////
  /// ADC/DAC for each ADC range 
  vector<double> inlADCPerDAC;          
  /// sigma/val for all vals in INT_NONLIN group
  float inlSigPct;             

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
