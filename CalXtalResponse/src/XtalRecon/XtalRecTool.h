#ifndef XtalRecTool_h
#define XtalRecTool_h
// $Header$

// LOCAL
#include "CalXtalResponse/IXtalRecTool.h"

// GLAST
#include "CalUtil/CalDefs.h"
#include "CalUtil/CalArray.h"

// EXTLIB
#include "GaudiKernel/AlgTool.h"

// STD
#include <cstring>

/*! @class XtalRecTool
  \author Zachary Fewtrell
  \brief Simple implementation of IXtalRecTool.  Faithfully pasted
  from XtalRecAlg v5r6p1

*/

class Event::CalDigi;
class ICalCalibSvc;
class IGlastDetSvc;
class Point;
class TTree;
class TFile;

class XtalRecTool : public AlgTool, 
                    virtual public IXtalRecTool 
{
 public:

  /// default ctor, declares jobOptions.
  XtalRecTool( const string& type, 
               const string& name, 
               const IInterface* parent);

  /// gets needed parameters and pointers to required services
  StatusCode initialize();

  StatusCode finalize();

  StatusCode calculate(const Event::CalDigi &digi,
                       Event::CalXtalRecData &xtalRec,
					   CalUtil::CalArray<CalUtil::FaceNum, bool> &belowThresh,
                       bool &xtalBelowThresh,
                       CalUtil::CalArray<CalUtil::FaceNum, bool> &saturated,
					   const INeighborXtalkTool *xtalTool=0,
                       const Event::EventHeader *evtHdr=0
                       );


 private:
  /** \brief convert large diode CIDAC scale to small diode CIDAC scale
      for given xtal face & pos 

      
      \param face which xtal face to process (uses current xtal from
      m_dat struct)

      \param pos longitudinal position along xtal (mm)
      \param largeCIDAC input large diode CIDAC value
      \param smallCIDAC output small diode CIDAC value
  */
  StatusCode largeCIDAC2Small(CalUtil::FaceNum face, float pos, float largeCIDAC, 
                            float &smallCIDAC);


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
  auto_ptr<TFile> m_tupleFile;


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
	CalUtil::CalArray<CalUtil::FaceNum, unsigned short> adc;

    CalUtil::CalArray<CalUtil::FaceNum, float> adcPed;

    float  ene;
    CalUtil::CalArray<CalUtil::FaceNum, float> faceSignal;
    CalUtil::CalArray<CalUtil::DiodeNum, float> asymCtr;

    float  pos;
    CalUtil::CalArray<CalUtil::FaceNum, float> cidac;

    float  asym;
    float  meanCIDAC;

    /// calibration constant
    CalUtil::CalArray<CalUtil::FaceNum, float> ped;
    /// calibration constant
    CalUtil::CalArray<CalUtil::FaceNum, float> pedSig;
    /// calibration constant
    CalUtil::CalArray<CalUtil::FaceNum, float> lacThresh;
    /// calibration constant
    CalUtil::CalArray<CalUtil::DiodeNum, float> mpd;

    /// calibration constant
    CalUtil::CalArray<CalUtil::FaceNum, float> h1Limit;

    CalUtil::TwrNum  twr;
    CalUtil::LyrNum  lyr;
    CalUtil::ColNum  col;

    CalUtil::CalArray<CalUtil::FaceNum, CalUtil::RngNum> rng;

    CalUtil::CalArray<CalUtil::FaceNum, char> belowThresh;

    char  xtalBelowThresh;

    CalUtil::CalArray<CalUtil::FaceNum, char> saturated;

    CalUtil::CalArray<CalUtil::FaceNum, CalUtil::DiodeNum> diode;

    CalUtil::XtalIdx xtalIdx;

  };

  /// used for each tuple.Fill() operation
  AlgData m_dat;
};


#endif
