#ifndef XtalRecTool_h
#define XtalRecTool_h

// LOCAL
#include "CalXtalResponse/IXtalRecTool.h"
#include "CalXtalResponse/ICalCalibSvc.h"

// GLAST
#include "CalUtil/CalDefs.h"
#include "CalUtil/CalArray.h"
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"

// EXTLIB
#include "GaudiKernel/AlgTool.h"
#include "TTree.h"

// STD
#include <cstring>

using namespace CalUtil;
using namespace Event;

/*! @class XtalRecTool
  \author Zachary Fewtrell
  \brief Simple implementation of IXtalRecTool.  Faithfully pasted
  from XtalRecAlg v5r6p1

*/

class XtalRecTool : 
public AlgTool, 
virtual public IXtalRecTool 
{
 public:

  /// default ctor, declares jobOptions.
  XtalRecTool( const string& type, 
               const string& name, 
               const IInterface* parent);

  /// gets needed parameters and pointers to required services
  virtual StatusCode initialize();

  virtual StatusCode finalize();

  StatusCode calculate(const Event::CalDigi &digi,
                       Event::CalXtalRecData &xtalRec,
                       CalArray<FaceNum, bool> &belowThresh,
                       bool &xtalBelowThresh,
                       CalArray<FaceNum, bool> &saturated,
                       CalTupleEntry *calTupleEnt=0,
                       const Event::EventHeader *evtHdr=0
                       );


 private:
  /** \brief convert large diode DAC scale to small diode DAC scale
      for given xtal face & pos 

      
      \param face which xtal face to process (uses current xtal from
      m_dat struct)

      \param pos longitudinal position along xtal (mm)
      \param largeDAC input large diode DAC value
      \param smallDAC output small diode DAC value
  */
  StatusCode largeDAC2Small(FaceNum face, float pos, float largeDAC, 
                            float &smallDAC);


  /** \brief convert scalar position in mm from xtal center (longitudinal
      to 3d vector position
  */
  void pos2Point(float pos, Point &pXtal);

  /// retrieve needed calibration constants for this xtal load into m_dat
  StatusCode retrieveCalib();

  /// name of CalCalibSvc to use for calib constants.
  StringProperty m_calCalibSvcName;                         
  /// pointer to CalCalibSvc object.
  ICalCalibSvc *m_calCalibSvc;  

  /// pointer to the Glast Detector Service
  IGlastDetSvc* m_detSvc; 

  /** \brief filename of XtalRecToolTuple. No file created if set to default=""
   */
  StringProperty m_tupleFilename;
  /// pointer to XtalRecToolTuple (TTree actually).  tuple is ignored
  /// if pointer is NULL
  TTree *m_tuple;
  /// pointer to XtalRecToolTuple file.
  TFile *m_tupleFile;


  /// length of CsI xtal in mm
  float m_CsILength;
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

    unsigned   RunID;
    unsigned   EventID;
    CalArray<FaceNum, unsigned short> adc;

    CalArray<FaceNum, float> adcPed;

    float  ene;
    CalArray<FaceNum, float> faceSignal;
    CalArray<DiodeNum, float> asymCtr;

    float  pos;
    CalArray<FaceNum, float> dac;

    float  asym;
    float  meanDAC;

    /// calibration constant
    CalArray<FaceNum, float> ped;
    /// calibration constant
    CalArray<FaceNum, float> pedSig;
    /// calibration constant
    CalArray<FaceNum, float> lacThresh;
    /// calibration constant
    CalArray<DiodeNum, float> mpd;

    /// calibration constant
    CalArray<FaceNum, float> h1Limit;


    TwrNum  twr;
    LyrNum  lyr;
    ColNum  col;

    CalArray<FaceNum, RngNum> rng;

    CalArray<FaceNum, char> belowThresh;

    char  xtalBelowThresh;

    CalArray<FaceNum, char> saturated;

    CalArray<FaceNum, DiodeNum> diode;

    XtalIdx xtalIdx;

  };

  /// used for each tuple.Fill() operation
  AlgData m_dat;
};


#endif
