#ifndef XtalRecTool_h
#define XtalRecToo_h

// LOCAL
#include "CalXtalResponse/IXtalRecTool.h"
#include "CalXtalResponse/ICalCalibSvc.h"

// GLAST
#include "CalUtil/CalDefs.h"
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"

// EXTLIB
#include "GaudiKernel/AlgTool.h"
#include "TTree.h"

// STD
#include <cstring>

using namespace CalDefs;
using namespace Event;

/*! @class XtalRecTool
  \author Zachary Fewtrell
  \brief Simple implementation of IXtalRecTool.  Faithfully pasted from XtalRecAlg v5r6p1
*/

class XtalRecTool : 
public AlgTool, 
virtual public IXtalRecTool 
{
 public:

  /// default ctor, declares jobOptions.
  XtalRecTool::XtalRecTool( const string& type, 
                            const string& name, 
                            const IInterface* parent);

  /// gets needed parameters and pointers to required services
  virtual StatusCode initialize();

  virtual StatusCode finalize();

  StatusCode calculate(const Event::EventHeader &evtHdr,
                       const CalDigi &digi,
                       CalXtalRecData &xtalRec,
                       bool &belowThreshP,  
                       bool &belowThreshN,  
                       bool &xtalBelowThresh,
                       bool &saturatedP,     
                       bool &saturatedN,     
                       CalTupleEntry *calTupleEnt
                       );

 private:
  StatusCode largeDAC2Small(FaceNum face, float largeDAC, float &smallDAC);

  /** \brief convert scalar position in mm from xtal center (longitudinal
      to 3d vector position
  */
  void pos2Point(float pos, Point &pXtal);

  /// retrieve needed calibration constants for this xtal load into m_dat
  StatusCode retrieveCalib(CalXtalId xtalId);


  /// name of CalCalibSvc to use for calib constants.
  StringProperty m_calCalibSvcName;                         
  /// pointer to CalCalibSvc object.
  ICalCalibSvc *m_calCalibSvc;  

  /// pointer to the Glast Detector Service
  IGlastDetSvc* m_detSvc; 

  /** \brief filename of XtalRecToolTuple.  No file created if set to default=""          
   */
  StringProperty m_tupleFilename;
  /// pointer to XtalRecToolTuple (TTree actually).  tuple is ignored if pointer is NULL
  TTree *m_tuple;
  /// pointer to XtalRecToolTuple file.
  TFile *m_tupleFile;
  /// pointer to main branch in XtalRecTuple
  TBranch *m_tupleBranch;

  /// length of CsI xtal in mm
  double m_CsILength;
  /// the value of fLATObjects field, defining LAT towers 
  int m_eLATTowers; 
  
  /// the value of fTowerObject field, defining calorimeter module 
  int m_eTowerCAL;
  /// the value of fCellCmp field defining CsI crystal
  int m_eXtal;      
  /// number of geometric segments per Xtal
  int m_nCsISeg; 

  /** \brief holds local vars for current iteration of algorithm

  also used to populate hitTuple
  */

  struct AlgData {
    void Clear() {memset(this,0,sizeof(AlgData));}

    UInt_t   RunID;
    UInt_t   EventID;
    UShort_t adc[2];
    Float_t  adcPed[2]; 
    Float_t  ene;
    Float_t  faceSignal[2];
    Float_t  asymCtr[2];
    Float_t  pos;
    Float_t  dac[2];
    Float_t  asym;
    Float_t  meanDAC;

    /// calibration constant
    Float_t  ped[2];
    /// calibration constant
    Float_t  pedSig[2];
    /// calibration constant
    Float_t  lacThresh[2];
    /// calibration constant
    Float_t  mpd[2];
    /// calibration constant
    Float_t  h1Limit[2];

    // Move char types to end of structure to help ROOT treealignment issues.
    UChar_t  twr;
    UChar_t  lyr;
    UChar_t  col;
    UChar_t  success;
    UChar_t  rng[2];
    UChar_t  belowThresh[2];
    UChar_t  xtalBelowThresh;
    UChar_t  saturated[2];
    UChar_t  diode[2];
  };

  /// used for each tuple.Fill() operation
  AlgData m_dat;

  /// official ROOT tuple format string
  static const char *m_tupleDesc;
};


#endif
