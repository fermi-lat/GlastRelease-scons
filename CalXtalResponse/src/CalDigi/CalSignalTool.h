#ifndef CalSignalTool_H
#define CalSignalTool_H
// $Header$
/** @file     
    @author Z.Fewtrell

*/

// LOCAL 

// GLAST 
#include "CalXtalResponse/ICalSignalTool.h"


// EXTLIB
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/IIncidentListener.h"

// STD
#include <string>
#include <vector>

class IDataProviderSvc;
class IXtalSignalTool;
class IPrecalcCalibTool;
class IGlastDetSvc;


namespace Event {
  class McIntegratingHit;
};
class ICalCalibSvc;


/** @class CalSignalTool
 * \author  Z.Fewtrell
 * 
 * 
 * - Convert McIntegratingHitCol into diode signal levels (CIDAC) for each crystal diode.
 * - Cache results for remaindder of event to reduce redundant processing.
 *
 * jobOptions
 *  - CalCalibSvc (default="CalCalibSvc") - calibration data source
 *  - XtalSignalToolName (default="XtalSignalTool") - convert each McIntegratingHit into crystal signal levels
 *  - enableNoise (default=true) - enable poissonic (n_electrons)  and electronic (pedestal width) noise simulations
 *  - PrecalcCalibTool (default="PrecalcCalibTool") - stores calibrations 'derived' from base calibrations
 *
 */

class CalSignalTool : 
  public AlgTool, 
  virtual public ICalSignalTool,
  virtual public IIncidentListener {
  
public:
  /// default ctor, declares jobOptions
  CalSignalTool( const std::string& type, 
                 const std::string& name, 
                 const IInterface* parent);
  
  StatusCode initialize();
  StatusCode finalize () {return StatusCode::SUCCESS;}

  /// get signal level for given cal diode.
  StatusCode getDiodeSignal(CalUtil::DiodeIdx diodeIdx, float &signal);

  /// get signal level for given cal diode.
  StatusCode getTrigDiodeSignal(CalUtil::DiodeIdx diodeIdx, float &signal);
 
  /// get map bewteen McIntegrating hits and crystals
  const CalRelationMap *getCalRelationMap();

  /// hook the BeginEvent so that we can check our validity once per event.
  void handle ( const Incident& inc );

private:
  
  /// map diodes to electronic signal levels
  typedef CalUtil::CalVec<CalUtil::DiodeIdx, float> CalSignalMap;

  /// clear internal tables
  void newEvent();

  /// ensure that internal data is either sync'd with current event or empty
  StatusCode syncData();

  /// load internal maps with MC data from TDS
  StatusCode loadSignalMaps();

  /// apply electronic noise calcuation to entire cal
  StatusCode calcNoise();

  /// apply poissonic noise (based on electron count) calculation to single crystal
  StatusCode calcPoissonicNoiseXtal(const CalUtil::XtalIdx xtalIdx, CalSignalMap& signalMap);

  /// apply electronic noise (pedestal width) calculation to single crystal
  StatusCode calcElectronicNoiseXtal(const CalUtil::XtalIdx xtalIdx);

  /// register McHit <> Xtal relation
  StatusCode registerHitRel(Event::McIntegratingHit &hit);

  /// retrieve constants from db
  StatusCode retrieveConstants();

  /// sum hit signal to all applicable diodes in private store
  StatusCode sumHit(const Event::McIntegratingHit &hit);

  /// sum hit signal for the trigger to all applicable diodes in private store
  StatusCode sumHit4Trig(const Event::McIntegratingHit &hit);

  /////////////////////////
  //-- PRIVATE MEMBERS --//
  /////////////////////////
  
  /// map crystals to a bit mask for determining energy source
  typedef CalUtil::CalVec<CalUtil::XtalIdx, unsigned short> CalSourceMap;

  /// enumerate some source possibilities
  enum CalSourceType { simulation = 1 << 0,
                       overlay    = 1 << 1
                     };

  /// map diodes to electronic signal levels
  CalSignalMap m_calSignalMap;

  /// map diodes to electronic signal levels for trigger
  CalSignalMap m_calTrigSignalMap;

  /// map crystals to McIntegratingHit source
  CalSourceMap m_calSourceMap;

  /// map crystals to MCIntegratingHits
  CalRelationMap m_calRelMap;

  /// if current private data store is valid
  bool m_isValid;

  /// ptr to event svc
  IDataProviderSvc* m_evtSvc;

  /// volume ID enumeration
  int m_eTowerCAL;  

  /// volume ID enumeration
  int m_eLATTowers;

  /// volume ID enum
  int m_eMeasureX;

  /// volume ID enum
  int m_eXtal;

  /// gain - electrons/MeV 1=Sm, 0=Large
  int m_ePerMeV[2];       

  /// parameter to control amount of electronic noise (WBA - 2010)
  float m_electronicNoiseGain;

  /// name of Tool for calculating single xtal signal response.
  StringProperty m_xtalSignalToolName;

  /// pointer to xtal signal tool
  IXtalSignalTool* m_xtalSignalTool;

  /// enable noise simulation
  BooleanProperty m_enableNoise;

  /// enable crystal noise simulation
  BooleanProperty m_enableXtalNoise;

  /// enable electronics noise simulation
  BooleanProperty m_enableElecNoise;

  /// name of precalc calib tool
  StringProperty m_precalcCalibName;
  
  /// pointer to precalcCalibTool
  IPrecalcCalibTool *m_precalcCalib;


  /// name of CalCalibSvc to use for calib constants.
  StringProperty  m_calCalibSvcName;      
  
  /// pointer to CalCalibSvc object.
  ICalCalibSvc   *m_calCalibSvc; 

  /// used for constants & conversion routines.
  IGlastDetSvc* m_detSvc;


  /// list of active tower bays, populated at run time
  std::vector<CalUtil::TwrNum> m_twrList;

};

#endif // CalSignalTool_H
