/*
@file TrgConfigSvc.h

@brief provides Information about the Configuration
@author Eric Charles
   From Martin Kocian's TrgConfigSvc

$Header$

*/
#ifndef TrgConfigSvc_H
#define TrgConfigSvc_H 

// Include files

#include "ConfigSvc/IConfigSvc.h"
#include "GaudiKernel/Service.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/IIncidentListener.h"
#include "GaudiKernel/Property.h"

#include "configData/db/TrgConfigDB.h"

#include "enums/Lsf.h"
#include <string>


class IMootSvc;
class TrgConfig;
class FswEfcSampler;

/** @class ConfigSvc
* @brief Service to retrieve the LAT and Filter configuration
*
* Author:  E. Charles, from M. Kocian's TrgConfigSvc
*
*/

class ConfigSvc : virtual public Service, virtual public IConfigSvc  {

protected:
  
  static unsigned filterSlotKey(unsigned lpaMode, unsigned handlerId) {
    return (lpaMode << 8) | handlerId;
  }

  enum ConfigPart {
    // First FSW filters, use the handler IDs
    GAMMA = 1,
    MIP = 3,
    HIP = 4,
    DGN = 5,
    // Now LATC, starting at 8;
    GEM = 8,
    GEM_DFT = 9,
    ROI = 10
  };

public:

  /// C'tor, give it a name and a way to find service it needs
  ConfigSvc(const std::string& name, ISvcLocator* pSvcLocator); 
  
  /// Get and link up with the services we need
  StatusCode initialize();
  
  /// Finalize just cleans up
  StatusCode finalize();
     
  /// queryInterface - for implementing a Service this is necessary
  StatusCode queryInterface(const InterfaceID& riid, void** ppvUnknown);
  
  /// return the service type
  const InterfaceID& type() const;

  /// Make sure we have the correct interfaceID
  static const InterfaceID& interfaceID() {
    return IConfigSvc::interfaceID(); 
  }

  /* 
   * Stuff from the IConfigSvc interface
   */
  
  /// get the MOOT key from MootSvc
  virtual unsigned getMootKey() const;

  /// get the GEM configuration object
  virtual const TrgConfig* getTrgConfig() const;

  /// get the information about the prescalers
  virtual const FswEfcSampler* getFSWPrescalerInfo( enums::Lsf::Mode mode, unsigned handlerId, 
						    unsigned int& fmxKey ) const;      

protected:

  /// Get the full path of a file
  void getFullPath( const std::string& mootPath, std::string& fullPath ) const;
 
  /// check to see if we have a new MOOT key
  bool newMootKey() const;
 
  /// clear our the read data strutures
  void clearCache() const;
  
  /// read the TrgConfig
  bool readTrgConfig( const std::string& gemFileName, 
		      const std::string& gemDftFileName, 
		      const std::string& roiFileName ) const;    

  /// function to get the prescale factors
  FswEfcSampler* getFswSampler(unsigned lpaMode, unsigned handlerId, unsigned int& fmxKey ) const;

  /// function to get the cached prescale factors
  FswEfcSampler* getFswSamplerFromCache(unsigned lpaMode, unsigned handlerId, unsigned int& fmxKey ) const;
  
  /// read an xml file with the prescaler info 
  FswEfcSampler* readEfcFromFile(unsigned mode, unsigned handlerId, 
				 const std::string& fileName, unsigned fmxKey) const;
  
private:
  
  StringProperty m_gammaFilterXmlFile;
  StringProperty m_dgnFilterXmlFile;
  StringProperty m_mipFilterXmlFile;
  StringProperty m_hipFilterXmlFile;    
  StringProperty m_gemXmlFile;
  StringProperty m_gemDftXmlFile;
  StringProperty m_roiXmlFile;
  
  IMootSvc*      m_mootSvc;

  mutable TrgConfig*                               m_trgConfig;
  mutable std::map<unsigned,FswEfcSampler*>        m_fswEfcCache;  
  mutable std::map<unsigned,unsigned>              m_fswCdmKeyCache;

  mutable unsigned                                 m_noMOOTMask;
  mutable unsigned                                 m_mootKey;
  mutable unsigned                                 m_latcKey;
};


#endif // TrgConfigSvc_H

