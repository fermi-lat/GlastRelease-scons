#ifndef XtalSignalTool_h
#define XtalSignalTool_h
//  $Header$

// LOCAL
#include "IXtalSignalTool.h"

// GLAST
#include "CalUtil/CalDefs.h"
#include "CalUtil/CalArray.h"


// EXTLIB
#include "GaudiKernel/AlgTool.h"
#include "TFile.h"

// STD
#include <cstring>

// FORWARD DECLARATIONS
class ICalCalibSvc;
class TTree;
class IDataProviderSvc;

/*! \class XtalSignalTool
  \author Zachary Fewtrell
  \brief Default implementation of IXtalSignalTool.  

  convert McIntegratingHit into diode signal for single xtal and hit

  jobOptions:
  - CalCalibSvc (default="CalCalibSvc") - calibration data source
  - tupleFile name (default="" - disabled) - generate debugging tuple file when string is non-zero length 
  
*/


class XtalSignalTool : public AlgTool, 
                       virtual public IXtalSignalTool
{
public:
  /// default ctor, declares jobOptions
  XtalSignalTool( const string& type, 
                  const string& name, 
                  const IInterface* parent);

  /// gets needed parameters and pointers to required services
  StatusCode initialize();

  StatusCode finalize();

  StatusCode calculate(const Event::McIntegratingHit &hitList,
                       CalUtil::CalArray<CalUtil::XtalDiode, float> &cidacArray);


private:
  /// take CsI MC integrated hit & sum deposited energy to all xtal diodes
  StatusCode sumCsIHit(const Event::McIntegratingHit &hit,
                       CalUtil::CalArray<CalUtil::XtalDiode, float> &cidacArray);

  /// calculate signal level from direct diode deposit 
  StatusCode sumDiodeHit(const Event::McIntegratingHit &hit,
                         CalUtil::CalArray<CalUtil::XtalDiode, float> &cidacArray);
  
  /// convert mcHit to longitudinal mm from xtal center
  float hit2pos(const Event::McIntegratingHit &hit);

  /// retrieve constants from Db
  StatusCode retrieveConstants();

  /// initialize optional tuple
  StatusCode initTuple();

  /// name of CalCalibSvc to use for calib constants.
  StringProperty  m_calCalibSvcName;      
  
  /// pointer to CalCalibSvc object.
  ICalCalibSvc   *m_calCalibSvc; 

  //-- Gaudi supplId constants --//
  /// number of geometric segments per Xtal
  int m_nCsISeg;                         

  /// Xtal length
  float m_CsILength;                    
  /// indicates volume inside CsI crystal
  int m_eXtal;
  /// detModel identifier for sm minus-side diode
  int m_eDiodeMSm;                       
  /// detModel identifier for sm plus-side diode
  int m_eDiodePSm;                       
  /// detModel identifier for large minus-side diode
  int m_eDiodeMLarge;                    
  /// detModel identifier for large plus-side diode
  int m_eDiodePLarge;                    
  
  /// used for calculating direct diode deposits.
  static const float m_ePerMeVInDiode;
  /// gain - electrons/MeV 1=Sm, 0=Large
  int m_ePerMeV[2];                      

  /// gain for diode hits
  
  /// filename of XtalSignalToolTuple.  No file created if set to default=""
  StringProperty m_tupleFilename;
  /// pointer to XtalSignalToolTuple (TTree actually).  tuple is ignored
  /// if pointer is NULL
  TTree *m_tuple;
  /// pointer to XtalSignalToolTuple file.
  auto_ptr<TFile> m_tupleFile;

  /** \brief holds local vars for current iteration of algorithm

  also used to populate debuggin tuple
  */
  struct AlgData {
    void Clear() {memset(this,0,sizeof(AlgData));}

    unsigned   RunID;
    unsigned   EventID;

    /// cidac values for each adc range
    CalUtil::CalArray<CalUtil::XtalDiode, float> diodeCIDAC;

    /// current xtal
    CalUtil::XtalIdx xtalIdx;

    /// calibration constant
    CalUtil::CalArray<CalUtil::DiodeNum, float> mpd;

  };

  AlgData m_dat;

  /// ptr to event svc
  IDataProviderSvc* m_evtSvc;



};


#endif
