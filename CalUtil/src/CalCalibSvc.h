#ifndef CalCalibSvc_H
#define CalCalibSvc_H 1

// Include files
// ANSI C++
#include <sstream>

// GAUDI
#include "GaudiKernel/Service.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/IIncidentListener.h"

// ROOT
#include "THashList.h"

// GLAST
#include "CalibData/Cal/CalCalibGain.h"
#include "CalibData/Cal/CalCalibIntNonlin.h"
#include "CalibData/Cal/CalCalibLightAsym.h"
#include "CalibData/Cal/CalCalibLightAtt.h"
#include "CalibData/Cal/CalCalibMuSlope.h"
#include "CalibData/Cal/CalCalibPed.h"

#include "CalUtil/ICalCalibSvc.h"


/** @class CalCalibSvc
 * \brief Instatiates ICalCalibSvc interface, retrieves data from CalibDataSvc
 *
 * handles:
 * - data storage/destruction
 * - communication with Gleam lower level services
 * - checking of data validity period  
 * - extraction of cal-specific constants out of generic data objects
 * - creation/caching of spline function objects where needed.
 *
 * \author  Zachary Fewtrell
 *
 */

class CalCalibSvc : public Service, virtual public ICalCalibSvc, virtual public IIncidentListener {
    
public:

  CalCalibSvc(const std::string& name, ISvcLocator* pSvcLocator); 

  StatusCode initialize();
  StatusCode execute();
  StatusCode finalize();

  /// queryInterface - for implementing a Service this is necessary
  StatusCode queryInterface(const IID& riid, void** ppvUnknown);

  static const InterfaceID& interfaceID() {
	 return ICalCalibSvc::interfaceID(); 
  }

  /// return the service type
  const IID& type() const;

  /// retrieve specified CAL_ElecGain value from CalibDataSvc
  StatusCode getGain(const idents::CalXtalId &xtalId, 
							idents::CalXtalId::XtalFace face,
							idents::CalXtalId::AdcRange range, 
							float &gain,
							float &sig);

  /// retrieve specified CAL_IntNonlin values from CalibDataSvc
  StatusCode getIntNonlin(const idents::CalXtalId &xtalId, 
								  idents::CalXtalId::XtalFace face,
								  idents::CalXtalId::AdcRange range, 
								  const std::vector< float > *&vals,
								  const std::vector< unsigned > *&dacs,
                          float &error);

  /// retrieve specified spline based on CAL_IntNonlin values form CalibDataSvc
  StatusCode getIntNonlin(const idents::CalXtalId &xtalId,
                          idents::CalXtalId::XtalFace face,
                          idents::CalXtalId::AdcRange range,
                          const TSpline3 *&intNonlinSpline);

  
  /// retrieve specified CAL_LightAsymetry values form CalibDataSvc
  StatusCode getLightAsym(const idents::CalXtalId &xtalId, 
								  idents::CalXtalId::XtalFace face,
								  idents::CalXtalId::AdcRange range, 
								  const std::vector< float > *&vals,
								  float &error);
  
  /// retrieve specified CAL_LightAtt values form CalibDataSvc
  StatusCode getLightAtt(const idents::CalXtalId &xtalId, 
								 idents::CalXtalId::XtalFace face,
								 idents::CalXtalId::AdcRange range, 
								 float &att,
								 float &norm);

  /// retrieve specified CAL_MuSlope values form CalibDataSvc
  StatusCode getMuSlope(const idents::CalXtalId &xtalId, 
								idents::CalXtalId::XtalFace face,
								idents::CalXtalId::AdcRange range, 
								float &slope,
								float &error);
  
  /// \brief retrieve specified CAL_Ped values form CalibDataSvc
  StatusCode getPed(const idents::CalXtalId &xtalId, 
						  idents::CalXtalId::XtalFace face,
						  idents::CalXtalId::AdcRange range, 
						  float &avr,
						  float &sig,
						  float &cosAngle);

private:
  ////////////////////////////////////////////////
  ////// PARAMETER CACHE MANAGEMENT //////////////
  ////////////////////////////////////////////////

  /// generate string name for indexing cache of TSpline3 objects 
  std::string &generateSplineName(std::string &grName,
                                 const std::string grType,
                                 const idents::CalXtalId &xtalId,
                                 idents::CalXtalId::XtalFace face,
                                 idents::CalXtalId::AdcRange range);

  StatusCode initIntNonlinCache();          ///< get fresh serial # for int nonlin spline cache
  StatusCode checkIntNonlinCache();         ///< check to see that intnonlin spline cache
  StatusCode clearIntNonlinCache();         ///< clear intNonlin spline cache 

  /// load specified TSpline3 object into IntNonlinCache
  StatusCode loadIntNonlinSpline(const idents::CalXtalId &xtalId,
                                idents::CalXtalId::XtalFace face,
                                idents::CalXtalId::AdcRange range);
  /// retrieve one spline from IntNonlinCache, load it if necessary
  StatusCode retrieveIntNonlinSpline(const idents::CalXtalId &xtalId,
                                    idents::CalXtalId::XtalFace face,
                                    idents::CalXtalId::AdcRange range,
                                    const TSpline3 *&pSpline);

    
  /// hook the BeginEvent so that we can check our validity once per event.
  void handle ( const Incident& inc );

  THashList  m_intNonlinList;            ///< ROOT data containers
  // used to compare against a new event's ser no to see if data has been updated
  int m_intNonlinSerNo;

  /////////////////////////////////////////////////
  ////// END PARAMETER CACHE MANAGEMENT ///////////
  /////////////////////////////////////////////////

  // JobOptions PROPERTIES
  StringProperty m_calibDataSvcName;     ///< name of CalibDataSvc, main data source

  StringProperty m_defaultFlavor;        ///<  default flavor for all calib types, unless otherwise specified.

  // calib_item specific flavors override defaultFlavor
  StringProperty m_flavorGain;           ///< calib flavor override for gain constants
  StringProperty m_flavorIntNonlin;      ///< calib flavor override for int-nonlin constants
  StringProperty m_flavorLightAsym;      ///< calib flavor override for light asymetry constants
  StringProperty m_flavorLightAtt;       ///< calib flavor override for light attenuation constants
  StringProperty m_flavorMuSlope;        ///< calib flavor override for muon slope constants
  StringProperty m_flavorPed;            ///< calib flavor override for pedestalconstants

  // GAUDI RESOURCES
  IService         *m_pCalibDataSvc;     ///< pointer to CalibDataSvc
  IDataProviderSvc *m_pDataProviderSvc;  ///< pointer to IDataProviderSvc interface of CalibDataSvc
  	 
  // TDS paths
  std::string m_elecGainPath;            ///< TCDS pathname for gain data 
  std::string m_intNonlinPath;           ///< TCDS pathname for intnonlin data
  std::string m_lightAsymPath;           ///< TCDS pathname for lightasym data
  std::string m_lightAttPath;            ///< TCDS pathname for light attenuation data
  std::string m_muSlopePath;             ///< TCDS pathname for muon slope data
  std::string m_pedPath;                 ///< TCDS pathname for pedestal data
};


#endif // CalCalibSvc_H

