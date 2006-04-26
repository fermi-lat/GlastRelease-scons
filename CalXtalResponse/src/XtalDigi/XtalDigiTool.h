#ifndef XtalDigiTool_h
#define XtalDigiTool_h
//  $Header$

// LOCAL
#include "CalXtalResponse/IXtalDigiTool.h"
#include "CalXtalResponse/ICalCalibSvc.h"
#include "CalXtalResponse/ICalTrigTool.h"
#include "../CalCalib/IPrecalcCalibTool.h"

// GLAST
#include "CalUtil/CalDefs.h"


// EXTLIB
#include "GaudiKernel/AlgTool.h"
#include "TFile.h"
#include "TTree.h"

// STD
#include <cstring>

/*! \class XtalDigiTool
  \author Zachary Fewtrell
  \brief Official implementation of IXtalDigiTool.  

  sums IntMcHits into digital response for one Cal xtalIdx.
  
*/

class XtalDigiTool : public AlgTool, 
             virtual public IXtalDigiTool
{
 public:
  /// default ctor, declares jobOptions
  XtalDigiTool( const string& type, 
                const string& name, 
                const IInterface* parent);

  /// gets needed parameters and pointers to required services
  StatusCode initialize();

  StatusCode finalize();

  /// calculate xtal Adc response for all rngs basedon collection of
  /// McHits in xtal & diode regions
  
  StatusCode calculate(const vector<const Event::McIntegratingHit*> &hitList,
                       const Event::EventHeader *evtHdr,            
                       Event::CalDigi &calDigi,
                       CalArray<FaceNum, bool> &lacBits,
                       CalArray<XtalDiode, bool> &trigBits,
                       Event::GltDigi *glt,
                       bool zeroSuppress);
 private:
  /// take CsI MC integrated hit & sum deposited energy to all xtal diodes
  StatusCode sumCsIHit(const Event::McIntegratingHit &hit);

  /// sum direct diode deposit energy
  StatusCode sumDiodeHit(const Event::McIntegratingHit &hit);
  
  /// select best adc range for both faces
  StatusCode rangeSelect();

  /// populate Digi TDS class
  StatusCode fillDigi(Event::CalDigi &calDigi);

  /// convert mcHit to longitudinal mm from xtal center
  float hit2pos(const Event::McIntegratingHit &hit);

  /// name of CalCalibSvc to use for calib constants.
  StringProperty  m_calCalibSvcName;      
  
  /// skip randomly generated noise & signal fluctuations, makes
  /// certain debug & testing comparisons easier.
  BooleanProperty m_noRandNoise;        
    
  /// pointer to CalCalibSvc object.
  ICalCalibSvc   *m_calCalibSvc; 

  /// name of ICalTrigTool member
  StringProperty m_calTrigToolName;

  /// ICalTrigTool for calculating trigger response.
  ICalTrigTool *m_calTrigTool;

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
  /// If true, only output tuple row if LAC=true (saves much disk
  /// space, default = true)
  BooleanProperty m_tupleLACOnly;
  /// pointer to XtalDigiToolTuple (TTree actually).  tuple is ignored
  /// if pointer is NULL
  TTree *m_tuple;
  /// pointer to XtalDigiToolTuple file.
  TFile *m_tupleFile;


  BooleanProperty m_allowNoiseOnlyTrig;

  /** \brief holds local vars for current iteration of algorithm

  also used to populate hitTuple
  */
  struct AlgData {
    void Clear() {memset(this,0,sizeof(AlgData));}

    unsigned   RunID;
    unsigned   EventID;

    /// ped subtracted adc.
    CalArray<XtalRng, float> adcPed;
    
    /// total # of MC int hits
    unsigned   nMCHits;
    /// # of MC int hits in CsI xtal
    unsigned   nCsIHits;
    /// # of MC int hits each diode
    CalArray<XtalDiode, unsigned> nDiodeHits;

    /// total ene deposited in xtal
    float  sumEneCsI;
    /// total ene from xtal and diodes
    float  sumEne;
    /// cidac values for each adc range
    CalArray<XtalDiode, float> diodeCIDAC;
    /// average energy deposit position weighted by ene of each deposit
    float  csiWeightedPos;
    
    /// calibration constant
    CalArray<XtalRng, float> ped;
    
    /// calibration constant CIDAC
    CalArray<XtalRng, float> pedSigCIDAC;

    /// trigger threholds in CIDAC units
    CalArray<XtalDiode, float> trigThreshCIDAC;

    /// calibration constant
    CalArray<FaceNum, float> lacThreshCIDAC;

    /// calibration constant
    CalArray<XtalRng, float> uldTholdADC;

    /// calibration constant
    CalArray<DiodeNum, float> mpd;

    CalArray<FaceNum, RngNum> rng;

    /// lac threshold flag
    CalArray<FaceNum, unsigned char> lac;
    
    CalArray<XtalDiode, unsigned char> trigBits;

    /// HEX1 range in xtal is saturated
    CalArray<FaceNum, unsigned char> saturated;

    XtalIdx xtalIdx;
  };
  AlgData m_dat;


  /// name of precalc calib tool
  StringProperty m_precalcCalibName;
  
  /// pointer to precalcCalibTool
  IPrecalcCalibTool *m_precalcCalib;

};


#endif
