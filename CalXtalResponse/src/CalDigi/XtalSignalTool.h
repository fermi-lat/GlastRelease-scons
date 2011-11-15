#ifndef XtalSignalTool_h
#define XtalSignalTool_h
//  $Header$
/** @file 
    @author Z.Fewtrell
*/

// LOCAL
#include "IXtalSignalTool.h"

// GLAST
#include "CalUtil/CalDefs.h"
#include "CalUtil/CalVec.h"

// EXTLIB
#include "GaudiKernel/AlgTool.h"

// STD
#include <cstring>

// FORWARD DECLARATIONS
class ICalCalibSvc;
class TTree;
class IDataProviderSvc;

/*! \class XtalSignalTool
  \author Z.Fewtrell
  \brief Default implementation of IXtalSignalTool.  

  convert McIntegratingHit into diode signal for single xtal and hit

  jobOptions:
  - CalCalibSvc (default="CalCalibSvc") - calibration data source
  
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

  /// calculate signal level on all diodes in single crystal from a
  /// group of McIntegratingHit energy deposits.
  StatusCode calculate(const Event::McIntegratingHit &hitList,
                       CalUtil::CalVec<CalUtil::XtalDiode, float> &cidacArray);


private:
  /// take CsI MC integrated hit & sum deposited energy to all xtal diodes
  StatusCode sumCsIHit(const Event::McIntegratingHit &hit,
                       CalUtil::CalVec<CalUtil::XtalDiode, float> &cidacArray);

  /// calculate signal level from direct diode deposit 
  StatusCode sumDiodeHit(const Event::McIntegratingHit &hit,
                         CalUtil::CalVec<CalUtil::XtalDiode, float> &cidacArray);
  
  /// convert mcHit to longitudinal mm from xtal center
  float hit2pos(const Event::McIntegratingHit &hit);

  /// retrieve constants from Db
  StatusCode retrieveConstants();

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
  
  /** \brief holds local vars for current iteration of algorithm

  also used to populate debuggin tuple
  */
  struct AlgData {
    void Clear() {
		fill(diodeCIDAC.begin(), diodeCIDAC.end(),0);
		xtalIdx = CalUtil::XtalIdx(); //< 0
		fill(mpd.begin(), mpd.end(),0);
	}

    /// cidac values for each adc range 
    CalUtil::CalVec<CalUtil::XtalDiode, float> diodeCIDAC;

    /// current xtal
    CalUtil::XtalIdx xtalIdx;

    /// calibration constant
    CalUtil::CalVec<CalUtil::DiodeNum, float> mpd;

  };

  AlgData m_dat;

};


#endif
