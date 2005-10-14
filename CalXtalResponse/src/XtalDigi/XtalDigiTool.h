#ifndef XtalDigiTool_h
#define XtalDigiTool_h

// LOCAL
#include "CalXtalResponse/IXtalDigiTool.h"
#include "CalXtalResponse/ICalCalibSvc.h"

// GLAST
#include "CalUtil/CalDefs.h"
#include "CalUtil/ICalFailureModeSvc.h"


// EXTLIB
#include "GaudiKernel/AlgTool.h"
#include "TTree.h"

// STD
#include <cstring>

/*! \class XtalDigiTool
  \author Zachary Fewtrell
  \brief Official implementation of IXtalDigiTool.  sums IntMcHits into digital response for one Cal xtalId.

*/

class XtalDigiTool : public AlgTool, virtual public IXtalDigiTool {
 public:
  /// default ctor, declares jobOptions
  XtalDigiTool::XtalDigiTool( const string& type, 
                              const string& name, 
                              const IInterface* parent);

  /// gets needed parameters and pointers to required services
  virtual StatusCode initialize();

  virtual StatusCode finalize();

  /// calculate xtal Adc response for all rngs basedon collection of McHits in xtal & diode regions
  StatusCode calculate(CalXtalId xtalId, 
                       const vector<const Event::McIntegratingHit*> &hitList,
                       const Event::EventHeader &evtHdr,            
                       Event::CalDigi &calDigi,     // output 
                       bool &lacP,                  // output 
                       bool &lacN,                  // output 
                       bool &fleP,                  // output 
                       bool &fleN,                  // output 
                       bool &fheP,                  // output 
                       bool &fheN                   // output 
                       );
 private:
  /// retrieve needed calibration constants for this xtal load into m_dat
  StatusCode retrieveCalib(CalXtalId xtalId);

  /// take CsI MC integrated hit & sum deposited energy to all xtal diodes
  StatusCode XtalDigiTool::sumCsIHit(CalXtalId xtalId, const Event::McIntegratingHit &hit);

  /// sum direct diode deposit energy
  StatusCode XtalDigiTool::sumDiodeHit(CalXtalId xtalId, const Event::McIntegratingHit &hit);
  
  /// simulate FLE & FHE triggers
  StatusCode XtalDigiTool::simTriggers();

  /// select best adc range for both faces
  StatusCode XtalDigiTool::rangeSelect();

  /// populate Digi TDS class
  StatusCode XtalDigiTool::fillDigi(CalXtalId xtalId, Event::CalDigi &calDigi);

  /// convert mcHit to longitudinal mm from xtal center
  float XtalDigiTool::hit2pos(const Event::McIntegratingHit &hit);

  /// name of CalCalibSvc to use for calib constants.
  StringProperty  m_calCalibSvcName;      
  
  /// \brief
  ///
  /// skip randomly generated noise & signal fluctuations, makes certain debug & testing comparisons easier.
  BooleanProperty m_noRandNoise;        
    
  /// pointer to CalCalibSvc object.
  ICalCalibSvc   *m_calCalibSvc; 

  /// pointer to failure mode service
  ICalFailureModeSvc* m_FailSvc;

  //-- Gaudi supplId constants --//
  /// number of geometric segments per Xtal
  int m_nCsISeg;                         
  int m_eXtal;
  /// Xtal length
  float m_CsILength;                    
  /// detModel identifier for sm minus-side diode
  int m_eDiodeMSm;                       
  /// detModel identifier for sm plus-side diode
  int m_eDiodePSm;                       
  /// detModel identifier for large minus-side diode
  int m_eDiodeMLarge;                    
  /// detModel identifier for large plus-side diode
  int m_eDiodePLarge;                    
  
  float m_ePerMeVInDiode;
  /// gain - electrons/MeV 1=Sm, 0=Large
  int m_ePerMeV[2];                      
  /// max val for ADC
  int m_maxAdc;  

  
  /// filename of XtalDigiToolTuple.  No file created if set to default=""
  StringProperty m_tupleFilename;
  /// If true, only output tuple row if LAC=true (saves much disk space, default = true)
  BooleanProperty m_tupleLACOnly;
  /// pointer to XtalDigiToolTuple (TTree actually).  tuple is ignored if pointer is NULL
  TTree *m_tuple;
  /// pointer to XtalDigiToolTuple file.
  TFile *m_tupleFile;
  /// pointer to main branch in XtalRecTuple
  TBranch *m_tupleBranch;

  BooleanProperty m_allowNoiseOnlyTrig;

  /** \brief holds local vars for current iteration of algorithm

  also used to populate hitTuple
  */
  struct AlgData {
    void Clear() {memset(this,0,sizeof(AlgData));}

    UInt_t   RunID;
    UInt_t   EventID;
    Float_t  adc[2][4];
	/// ped subtracted adc.
	Float_t  adcPed[2][4];
    
	/// total # of MC int hits
    UInt_t   nMCHits;
    /// # of MC int hits in CsI xtal
    UInt_t   nCsIHits;
    /// # of MC int hits each diode
    UInt_t   nDiodeHits[2][2];
    /// total ene deposited in xtal
    Float_t  sumEneCsI;
    /// total ene from xtal and diodes
    Float_t  sumEne;
    /// dac values for eachdiode
    Float_t  diodeDAC[2][2];
    /// average energy deposit position weighted by ene of each deposit
    Float_t  csiWeightedPos;
    
    /// calibration constant
    Float_t  ped[2][4];
    /// calibration constant
    Float_t  pedSig[2][4];
    /// calibration constant
    Float_t  lacThresh[2];
    /// calibration constant
    Float_t  fleThresh[2];
    /// calibration constant
    Float_t  fheThresh[2];
    /// calibration constant
    Float_t  uldThold[2][4];
    /// calibration constant
    Float_t  mpd[2];


	// Move char types to end of structure to help ROOT treealignment issues.
	UChar_t  twr;
    UChar_t  lyr;
    UChar_t  col;
    UChar_t  success;
    UChar_t  rng[2];
	/// lac threshold flag
    UChar_t   lac[2];
    /// fle trigger flag
    UChar_t   fle[2];
    /// fhe trigger flag
    UChar_t   fhe[2];
    /// HEX1 range in xtal is saturated
    UChar_t   saturated[2];
    
    
  };

  AlgData m_dat;
  
  static const char *m_tupleDesc;
};


#endif
