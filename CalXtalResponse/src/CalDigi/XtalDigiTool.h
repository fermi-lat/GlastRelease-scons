#ifndef XtalDigiTool_h
#define XtalDigiTool_h
//  $Header$

// LOCAL
#include "CalXtalResponse/IXtalDigiTool.h"

// GLAST
#include "CalUtil/CalDefs.h"
#include "CalUtil/CalArray.h"


// EXTLIB
#include "GaudiKernel/AlgTool.h"
#include "TFile.h"

// STD
#include <cstring> /// memset

// forward declarations
class ICalCalibSvc;
class TTree;
class ICalFailureModeSvc;
class IDataProviderSvc;
class IPrecalcCalibTool;

/*! \class XtalDigiTool
  \author Zachary Fewtrell
  \brief Default implementation of IXtalDigiTool.  

  Convert Diode signals to CalDigi object for one Cal crystal.

  jobOptions:
  - CalCalibSvc (default="CalCalibSvc") - Cal calibration data source
  - PrecalcCalibTool (deefault="PrecalcCalibTool") - Derived Cal calibration data source.
  - tupleFilename (default="" (disabled)) - optional debugging tuple file (leave set to "" for disabled)
  - tupleLACOnly (default=true) - disable tuple output for crystals below LAC threshold
  
  
*/

class XtalDigiTool : public AlgTool, 
             virtual public IXtalDigiTool
{
 public:
  /// default ctor, declares jobOptions
  XtalDigiTool( const std::string& type, 
                const std::string& name, 
                const IInterface* parent);

  /// gets needed parameters and pointers to required services
  StatusCode initialize();

  StatusCode finalize();

  StatusCode calculate(const ICalSignalTool::XtalSignalMap &cidac,
                       Event::CalDigi &calDigi,
                       CalUtil::CalArray<CalUtil::FaceNum, bool> &lacBits,
                       bool zeroSuppress
                       );
 private:
  
  /// select best adc range for both faces
  StatusCode rangeSelect();

  /// get failureMode bits
  StatusCode getFailureStatus(const CalUtil::CalArray<CalUtil::FaceNum, bool> &lacBits);

  /// populate Digi TDS class
  StatusCode fillDigi(Event::CalDigi &calDigi);

  /// name of CalCalibSvc to use for calib constants.
  StringProperty  m_calCalibSvcName;      
  
  /// pointer to CalCalibSvc object.
  ICalCalibSvc   *m_calCalibSvc; 

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
  auto_ptr<TFile> m_tupleFile;

  /** \brief holds local vars for current iteration of algorithm

  also used to populate hitTuple
  */
  struct AlgData {
    void Clear() {memset(this,0,sizeof(AlgData));}

    unsigned   RunID;
    unsigned   EventID;

    /// ped subtracted adc.
    CalUtil::CalArray<CalUtil::XtalRng, float> adcPed;
    
    /// cidac values for each adc range
    CalUtil::CalArray<CalUtil::XtalDiode, float> diodeCIDAC;
    
    /// calibration constant
    CalUtil::CalArray<CalUtil::XtalRng, float> ped;
    
    /// calibration constant
    CalUtil::CalArray<CalUtil::FaceNum, float> lacThreshCIDAC;

    /// calibration constant
    CalUtil::CalArray<CalUtil::XtalRng, float> uldTholdADC;

    CalUtil::CalArray<CalUtil::FaceNum, CalUtil::RngNum> rng;

    /// lac threshold flag
    CalUtil::CalArray<CalUtil::FaceNum, unsigned char> lac;
    
    /// HEX1 range in xtal is saturated
    CalUtil::CalArray<CalUtil::FaceNum, unsigned char> saturated;

    CalUtil::XtalIdx xtalIdx;
    
    unsigned short failureStatus;
  };

  AlgData m_dat;

  /// pointer to optional CalFailureMode
  ICalFailureModeSvc *m_calFailureModeSvc;

  /// ptr to event svc
  IDataProviderSvc* m_evtSvc;

  /// name of precalc calib tool
  StringProperty m_precalcCalibName;
  
  /// pointer to precalcCalibTool
  IPrecalcCalibTool *m_precalcCalib;
};


#endif
