#ifndef CalCalibSvc_H
#define CalCalibSvc_H 1

// Include files
// ANSI C++
#include <strstream>

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
 * @brief Instatiates ICalCalibSvc interface
 *
 * Author:  Z. Fewtrell
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

  StatusCode getGain(const idents::CalXtalId &xtalId, 
							idents::CalXtalId::XtalFace face,
							idents::CalXtalId::AdcRange range, 
							float &gain,
							float &sig);

  StatusCode getIntNonlin(const idents::CalXtalId &xtalId, 
								  idents::CalXtalId::XtalFace face,
								  idents::CalXtalId::AdcRange range, 
								  const std::vector< float > *&vals,
								  const std::vector< unsigned > *&dacs,
                          float &error);

  StatusCode getIntNonlin(const idents::CalXtalId &xtalId,
                          idents::CalXtalId::XtalFace face,
                          idents::CalXtalId::AdcRange range,
                          const TSpline3 *&intNonlinSpline);

  
  StatusCode getLightAsym(const idents::CalXtalId &xtalId, 
								  idents::CalXtalId::XtalFace face,
								  idents::CalXtalId::AdcRange range, 
								  const std::vector< float > *&vals,
								  float &error);
  
  StatusCode getLightAtt(const idents::CalXtalId &xtalId, 
								 idents::CalXtalId::XtalFace face,
								 idents::CalXtalId::AdcRange range, 
								 float &att,
								 float &norm);

  StatusCode getMuSlope(const idents::CalXtalId &xtalId, 
								idents::CalXtalId::XtalFace face,
								idents::CalXtalId::AdcRange range, 
								float &slope,
								float &error);
  
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

  StatusCode initIntNonlinCache();          // get fresh CalCalibData */ser_no
  StatusCode checkIntNonlinCache();         // check to see that intnonlin data is still valid
  StatusCode clearIntNonlinCache();         // duh  

  /// load specified TSpline3 object into IntNonlinCache
  StatusCode loadIntNonlinSpline(const idents::CalXtalId &xtalId,
                                idents::CalXtalId::XtalFace face,
                                idents::CalXtalId::AdcRange range);
  /// retrieve one face/range spline into IntNonlinCache, load it if necessary
  StatusCode retrieveIntNonlinSpline(const idents::CalXtalId &xtalId,
                                    idents::CalXtalId::XtalFace face,
                                    idents::CalXtalId::AdcRange range,
                                    const TSpline3 *&pSpline);

    
  // hook the BeginEvent so that we can keep parameter cache up to date
  void handle ( const Incident& inc );

  // ROOT data containers
  THashList  m_intNonlinList;
  // used to compare against a new event's ser no to see if data has been updated
  int m_intNonlinSerNo;

  /////////////////////////////////////////////////
  ////// END PARAMETER CACHE MANAGEMENT ///////////
  /////////////////////////////////////////////////

  // JobOptions PROPERTIES
  StringProperty m_calibDataSvcName;

  StringProperty m_defaultFlavor; // default flavor for all calib types

  // calib_item specific flavors override defaultFlavor
  StringProperty m_flavorGain;
  StringProperty m_flavorIntNonlin;
  StringProperty m_flavorLightAsym;
  StringProperty m_flavorLightAtt;
  StringProperty m_flavorMuSlope;
  StringProperty m_flavorPed;

  // GAUDI RESOURCES
  IService         *m_pCalibDataSvc;
  IDataProviderSvc *m_pDataProviderSvc;
  	 
  // TDS paths
  std::string m_elecGainPath;
  std::string m_intNonlinPath;
  std::string m_lightAsymPath;
  std::string m_lightAttPath;
  std::string m_muSlopePath;
  std::string m_pedPath;
};


#endif // CalCalibSvc_H

