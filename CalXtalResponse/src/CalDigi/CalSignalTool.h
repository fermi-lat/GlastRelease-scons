#ifndef CalSignalTool_H
#define CalSignalTool_H
// $Header$

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
 * \author  Zachary Fewtrell
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
  StatusCode execute();
  StatusCode finalize () {return StatusCode::SUCCESS;}

  /// return crystal signal level for given crystal, refer to McHits as needed (once per event)
  const CalSignalMap *getCalSignalMap();

  /// return crystal signal level for given crystal, refer to McHits as needed (once per event)
  StatusCode getXtalSignalMap(CalUtil::XtalIdx xtalIdx,
                              XtalSignalMap &xtalSignalMap);


  const CalRelationMap *getCalRelationMap();


  /// hook the BeginEvent so that we can check our validity once per event.
  void handle ( const Incident& inc );


private:
  /// clear internal tables (used by test app, unneeded during normal operations
  void newEvent();

  /// mark internal data state invalid & clear any existing data
  void invalidate();

  /// ensure that internal data is either sync'd with current event or empty
  StatusCode syncData();

  /// load internal maps with MC data from TDS
  StatusCode loadSignalMaps();

  /// apply electronic noise calcuation to entire cal
  StatusCode calcNoise();

  /// apply poissonic noise (based on electron count) calculation to single crystal
  StatusCode calcPoissonicNoiseXtal(const CalUtil::XtalIdx xtalIdx);

  /// apply electronic noise (pedestal width) calculation to single crystal
  StatusCode calcElectronicNoiseXtal(const CalUtil::XtalIdx xtalIdx);

  /// register McHit <> Xtal relation
  StatusCode registerHitRel(Event::McIntegratingHit &hit);

  /// retrieve constants from db
  StatusCode retrieveConstants();

  /// sum hit signal to all applicable diodes in private store
  StatusCode sumHit(const Event::McIntegratingHit &hit);

  /////////////////////////
  //-- PRIVATE MEMBERS --//
  /////////////////////////

  /// map diodes to electronic signal levels
  CalSignalMap m_calSignalMap;

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

  /// name of Tool for calculating single xtal signal response.
  StringProperty m_xtalSignalToolName;

  /// pointer to xtal signal tool
  IXtalSignalTool* m_xtalSignalTool;

  /// enable noise simulation
  BooleanProperty m_enableNoise;

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
