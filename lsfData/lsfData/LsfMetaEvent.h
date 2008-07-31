#ifndef LSFDATA_METAEVENT_H
#define LSFDATA_METAEVENT_H 1

#include <iostream>
//#include <map>

#include "lsfData/LsfTime.h"
#include "lsfData/LsfRunInfo.h"
#include "lsfData/LsfDatagramInfo.h"
#include "lsfData/LsfGemScalers.h"
#include "lsfData/LsfConfiguration.h"
#include "lsfData/LsfKeys.h"
#include "lsfData/LpaHandler.h"

/** @class MetaEvent
*
* $Header$
*/

namespace lsfData {
  
  class MetaEvent {
    
  public:
    
    MetaEvent( const RunInfo& run, const DatagramInfo& datagram, 
	       const GemScalers& scalers,
	       const Time& time,
	       const Configuration& configuration,
	       const LsfKeys& keys)
      :m_run(run),m_datagram(datagram),
       m_scalers(scalers),
       m_time(time),
       m_config(configuration.clone()),
       m_type(configuration.type()),
       m_keys(keys.clone()),
       m_ktype(keys.type()),
       m_gamma(0), m_pass(0), m_mip(0), m_hip(0), m_dgn(0), m_lpaHandler(0),
       m_mootKey(LSF_INVALID_UINT),
       m_mootAlias(""), m_compressionLevel(LSF_UNDEFINED),
       m_compressedSize(LSF_UNDEFINED) {
   
    }

    MetaEvent()
      :m_config(0),
       m_type(enums::Lsf::NoRunType),
       m_keys(0),
       m_ktype(enums::Lsf::NoKeysType),
       m_gamma(0), m_pass(0), m_mip(0), m_hip(0), m_dgn(0), m_lpaHandler(0),
       m_mootKey(LSF_INVALID_UINT),
       m_mootAlias(""),m_compressionLevel(LSF_UNDEFINED),
       m_compressedSize(LSF_UNDEFINED) {
    }
    
    MetaEvent( const MetaEvent& other ) :
       m_run(other.run()),
       m_datagram(other.datagram()),
       m_scalers(other.scalers()),
       m_time(other.time()),
       m_config(0),
       m_type(enums::Lsf::NoRunType),
       m_compressionLevel(other.compressionLevel()),
       m_compressedSize(other.compressedSize()) {
      if ( other.configuration() != 0 ) {
	m_config = other.configuration()->clone();
	m_type = other.configuration()->type();
      }
      if ( other.keys() != 0 ) {
	m_keys  = other.keys()->clone();
	m_ktype = other.keys()->type();
      }
      if (other.gammaFilter()) m_gamma = new lsfData::GammaHandler(*(other.gammaFilter()));
      if (other.mipFilter()) m_mip = new lsfData::MipHandler(*(other.mipFilter()));
      if (other.hipFilter()) m_hip = new lsfData::HipHandler(*(other.hipFilter()));
      if (other.dgnFilter()) m_dgn = new lsfData::DgnHandler(*(other.dgnFilter()));
      if (other.passthruFilter()) m_pass = new lsfData::PassthruHandler(*(other.passthruFilter()));
      if (other.lpaHandler()) m_lpaHandler = new lsfData::LpaHandler(*(other.lpaHandler()));
      m_mootKey = other.m_mootKey;
      m_mootAlias = other.m_mootAlias;
    }
    
    virtual ~MetaEvent(){
      delete m_config;
      delete m_keys;

      if (m_gamma) {
          delete m_gamma;
          m_gamma = 0;
      }
      if (m_mip) {
          delete m_mip;
          m_mip = 0;
      }
      if (m_hip) {
          delete m_hip;
          m_hip = 0;
      }
      if (m_dgn) {
          delete m_dgn;
          m_dgn = 0;
      }
      if (m_pass) {
          delete m_pass;
          m_pass = 0;
      }
      if (m_lpaHandler) {
          delete m_lpaHandler;
          m_lpaHandler = 0;
      }

    }

    inline void clear() {
       if (m_config) {
            delete m_config;
            m_config = 0;
        }
       if (m_keys) {
	 delete m_keys;
	 m_keys = 0;
       }
        m_run.clear();
        m_datagram.clear();
        m_scalers.clear();
        m_time.clear();      
        m_type = enums::Lsf::NoRunType;
	m_ktype = enums::Lsf::NoKeysType;
      if (m_gamma) {
          delete m_gamma;
          m_gamma = 0;
      }
      if (m_mip) {
          delete m_mip;
          m_mip = 0;
      }
      if (m_hip) {
          delete m_hip;
          m_hip = 0;
      }
      if (m_dgn) {
          delete m_dgn;
          m_dgn = 0;
      }
      if (m_pass) {
          delete m_pass;
          m_pass = 0;
      }
      if (m_lpaHandler) {
          delete m_lpaHandler;
          m_lpaHandler = 0;
      }
      m_mootKey = LSF_INVALID_UINT;
      m_mootAlias = "";

      m_compressionLevel = LSF_UNDEFINED;
      m_compressedSize = LSF_UNDEFINED;
    }

    /// Information about the run this event is from
    inline const RunInfo& run() const { return m_run; };

    /// Information about the datagram this event came in
    inline const DatagramInfo& datagram() const { return m_datagram; }

    /// The extended context records
    inline const GemScalers& scalers() const { return m_scalers; }

    /// Information about the time markers associated with this event
    inline const Time& time() const { return m_time; } 

    /// Information about the configuration keys associated with this event
    inline const Configuration* configuration() const { return m_config; }

    /// Translated configuration file keys for this event
    inline const LsfKeys* keys() const { return m_keys; };

    inline const MipHandler* mipFilter() const {
        return m_mip; }
    inline const HipHandler* hipFilter() const {
        return m_hip; }
    inline const DgnHandler* dgnFilter() const {
        return m_dgn;  }
    inline const PassthruHandler* passthruFilter() const {
        return m_pass;  }
    inline const GammaHandler* gammaFilter() const {
        return m_gamma;  }
    inline const LpaHandler* lpaHandler() const {
        return m_lpaHandler;  }


    inline unsigned int mootKey() const { return m_mootKey; }

    inline const std::string& mootAlias() const { return m_mootAlias; }

    inline int compressionLevel() const { return m_compressionLevel; }
    inline int compressedSize() const { return m_compressedSize; }

    /// set everything at once
    inline void set(const RunInfo& run, const DatagramInfo& datagram, 
		    const GemScalers& scalers,
		    const Time& time,
		    const Configuration& configuration,
		    const LsfKeys& keys) {
      m_run = run;
      m_datagram = datagram;
      m_scalers = scalers;
      m_time = time;
      if(m_config) delete m_config;
      m_config = configuration.clone();
      m_type = configuration.type();
      if(m_keys) delete m_keys;
      m_keys = keys.clone();
      m_ktype = keys.type();
    }

    // set the individual data members
    inline void setRun( const RunInfo& val) { m_run = val; };
    inline void setRunDataTransferId(unsigned int id) { m_run.setDataTransferId(id); };
    inline void setDatagram( const DatagramInfo& val) { m_datagram = val; };
    inline void setScalers( const GemScalers& val) { m_scalers = val; };
    inline void setTime( const Time& val) { m_time = val; }; 
    inline void setConfiguration( const Configuration& configuration ) {
      if (m_config) delete m_config;
      m_config = configuration.clone();
      m_type = configuration.type();
    }
    inline void setKeys( const LsfKeys& keys ) {
      if (m_keys) delete m_keys;
      m_keys = keys.clone();
      m_ktype = keys.type();
    }

    inline void setMootKey( unsigned int mootKey ) {
        m_mootKey = mootKey;
    }
  
    inline void setMootAlias( const char* mootAlias) {
        m_mootAlias = mootAlias;
    }

    inline void setCompressionLevel(int level) { m_compressionLevel=level;}
    inline void setCompressedSize(int size) { m_compressedSize = size; }

void addGammaHandler(const GammaHandler& gamma) {
    m_gamma = new lsfData::GammaHandler(gamma);
}
void addDgnHandler(const DgnHandler& dgn) {
    m_dgn = new lsfData::DgnHandler(dgn);
}
void addPassthruHandler(const PassthruHandler& pass) {
    m_pass = new lsfData::PassthruHandler(pass);
}
void addMipHandler(const MipHandler& mip) {
    m_mip = new lsfData::MipHandler(mip);
}
void addHipHandler(const HipHandler& hip) {
    m_hip = new lsfData::HipHandler(hip);
}
void addLpaHandler(const LpaHandler& lpa) {
    m_lpaHandler = new lsfData::LpaHandler(lpa);
}
    
  private:
    
    /// 
    RunInfo m_run;
    DatagramInfo m_datagram;
    GemScalers m_scalers;
    Time m_time;
    Configuration* m_config;
     
    enums::Lsf::RunType m_type;
    
    LsfKeys* m_keys;
    enums::Lsf::KeysType m_ktype;

    GammaHandler *m_gamma;   
    PassthruHandler *m_pass;
    MipHandler *m_mip;
    HipHandler *m_hip;
    DgnHandler *m_dgn;
    LpaHandler *m_lpaHandler;

    unsigned int m_mootKey;
    std::string  m_mootAlias;

    int m_compressionLevel;
    int m_compressedSize;

  };

}


#endif    // LSF_METAEVENT
