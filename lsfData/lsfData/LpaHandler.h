/** 
 * @class lsfData::LpaHandler
 *
 * @brief Class encapsulating event-handler information
 *
 * @author Heather Kelly <heather@slac.stanford.edu>
 *
 * $Header$
 */

#ifndef LSFDATA_LPAHANDLER_HH
#define LSFDATA_LPAHANDLER_HH


namespace lsfData {
    // forward declaration of sub-classes
 class LpaHandler;
 class DgnHandlerRsdV0;
 class GammaHandlerRsdV0;
 class HipHandlerRsdV0;
 class MipHandlerRsdV0;
 class PassthruHandlerRsdV0;


  /** class representing the result info from each active handler */
  class LpaHandler { 
  public:

    /** ctor and dump routine */
      LpaHandler():m_type(enums::Lsf::Unknown), m_masterKey(0xFFFFFFFF), m_cfgKey(0xFFFFFFFF), 
          m_cfgId(0xFFFFFFFF), m_state(enums::Lsf::INVALID), m_prescaler(enums::Lsf::UNSUPPORTED), m_version(0), 
          m_id(enums::Lsf::MaxHandlerIds), m_has(false) {};

      virtual ~LpaHandler() { };

      //LpaHandler* clone() const { std::cout << "LpaHandler::clone " << std::endl; return new LpaHandler(*this); };

    // override this to specify the HandlerId
    enums::Lsf::HandlerId id() const { return m_id; }    


    void print(const std::string &str="") const { }

    friend std::ostream& operator<< ( std::ostream& s, const LpaHandler& o) { return s; };
    
    /// This is a poor-man's dynamic cast
    //const LpaHandler* castToLpaHandler() const { return this; };

    void set(unsigned int masterKey, unsigned int cfgKey, unsigned int cfgId, 
            enums::Lsf::RsdState state, enums::Lsf::LeakedPrescaler prescaler, 
            unsigned int version, enums::Lsf::HandlerId id, bool has) {
            m_masterKey = masterKey;
            m_cfgKey = cfgKey;
            m_cfgId = cfgId;
            m_state = state;
            m_prescaler = prescaler;
            m_version = version;
            m_id = id;
            m_has = has;
        }
 
    unsigned int masterKey() const { return m_masterKey; };
    unsigned int cfgKey() const { return m_cfgKey; };
    unsigned int cfgId() const { return m_cfgId; };
    enums::Lsf::RsdState state() const { return m_state; };
    enums::Lsf::LeakedPrescaler prescaler() const { return m_prescaler; };
    unsigned int version() const { return m_version; };
    bool has() const { return m_has; }


      LpaHandler(const LpaHandler &other)  {
          m_type = other.m_type;
          m_masterKey = other.m_masterKey;
          m_cfgKey = other.m_cfgKey;
          m_cfgId = other.m_cfgId;
          m_state = other.m_state;
          m_prescaler = other.m_prescaler;
          m_version = other.m_version;
          m_id = other.m_id;
          m_has = other.m_has;
      } 

   // const char*                     typeName() const;
   // const char*                     handlerName() const;
  protected:
    enums::Lsf::HandlerType       m_type;
    unsigned int      m_masterKey;  /// FMX key of master CDM for handler (fixed for a run)
    unsigned int      m_cfgKey;     /// FMX key of current CDM for handler (may vary by mode)
    unsigned int      m_cfgId;      /// unique identifier of handler config (may vary by mode)
    enums::Lsf::RsdState   m_state;      /// overall filter result for the event
    enums::Lsf::LeakedPrescaler   m_prescaler;  /// Unsupported in RSD V0, otherwise indicates how event was leaked
    unsigned int     m_version;    /// Encoding version of handler-specific RSD
    enums::Lsf::HandlerId         m_id;         /// Handler identifier
    bool              m_has;        /// Indicates whether handler generated summary data


  };

  /** classes representing V0 of the handler-specific RSD's */
  class DgnHandlerRsdV0  {
  public:
      DgnHandlerRsdV0() {};
      ~DgnHandlerRsdV0() {};

    // DgnHandlerRsdV0* clone() const {
    //  std::cout << "Dgn clone" << std::endl;
    //  return new DgnHandlerRsdV0(*this);
   // }
    /// This is a poor-man's dynamic cast
   // const DgnHandlerRsdV0* castToDgnRsdV0() const { return this; };

    void setStatus(unsigned int status) { m_status = status; }
    unsigned int status() const { return m_status; }
    enums::Lsf::HandlerId id() const { return m_id; }    

    void set(unsigned int masterKey, unsigned int cfgKey, unsigned int cfgId, 
            enums::Lsf::RsdState state, enums::Lsf::LeakedPrescaler prescaler, 
            unsigned int version, enums::Lsf::HandlerId id, bool has) {
            m_masterKey = masterKey;
            m_cfgKey = cfgKey;
            m_cfgId = cfgId;
            m_state = state;
            m_prescaler = prescaler;
            m_version = version;
            m_id = id;
            m_has = has;
        }
 
    unsigned int masterKey() const { return m_masterKey; };
    unsigned int cfgKey() const { return m_cfgKey; };
    unsigned int cfgId() const { return m_cfgId; };
    enums::Lsf::RsdState state() const { return m_state; };
    enums::Lsf::LeakedPrescaler prescaler() const { return m_prescaler; };
    unsigned int version() const { return m_version; };
    bool has() const { return m_has; }


      DgnHandlerRsdV0(const DgnHandlerRsdV0 &other)  {
          m_type = other.m_type;
          m_masterKey = other.m_masterKey;
          m_cfgKey = other.m_cfgKey;
          m_cfgId = other.m_cfgId;
          m_state = other.m_state;
          m_prescaler = other.m_prescaler;
          m_version = other.m_version;
          m_id = other.m_id;
          m_has = other.m_has;
          m_status = other.m_status;
      } 

  private:
    enums::Lsf::HandlerType       m_type;
    unsigned int      m_masterKey;  /// FMX key of master CDM for handler (fixed for a run)
    unsigned int      m_cfgKey;     /// FMX key of current CDM for handler (may vary by mode)
    unsigned int      m_cfgId;      /// unique identifier of handler config (may vary by mode)
    enums::Lsf::RsdState   m_state;      /// overall filter result for the event
    enums::Lsf::LeakedPrescaler   m_prescaler;  /// Unsupported in RSD V0, otherwise indicates how event was leaked
    unsigned int     m_version;    /// Encoding version of handler-specific RSD
    enums::Lsf::HandlerId         m_id;         /// Handler identifier
    bool              m_has;        /// Indicates whether handler generated summary data
    unsigned int m_status;


  };
    
  class GammaHandlerRsdV0  {
  public:
      GammaHandlerRsdV0()  { };
      ~GammaHandlerRsdV0() { };

     //GammaHandlerRsdV0* clone() const {
     // std::cout << "Gamma clone" << std::endl;
     // return new GammaHandlerRsdV0(*this);
   // }
    /// This is a poor-man's dynamic cast
    //const GammaHandlerRsdV0* castToGammaRsdV0() const { return this; };


    void setStatus(unsigned int status, unsigned int stage, unsigned int energyValid, int energyInLeus) {
            m_status = status;
            m_stage = stage;
            m_energyValid = energyValid;
            m_energyInLeus = energyInLeus;
        }
 
    unsigned int status() const { return m_status; }
    unsigned int stage() const { return m_stage; }
    unsigned int energyValid() const { return m_energyValid; }
    int energyInLeus() const { return m_energyInLeus; }
    enums::Lsf::HandlerId id() const { return m_id; }    
 
    void set(unsigned int masterKey, unsigned int cfgKey, unsigned int cfgId, 
            enums::Lsf::RsdState state, enums::Lsf::LeakedPrescaler prescaler, 
            unsigned int version, enums::Lsf::HandlerId id, bool has) {
            m_masterKey = masterKey;
            m_cfgKey = cfgKey;
            m_cfgId = cfgId;
            m_state = state;
            m_prescaler = prescaler;
            m_version = version;
            m_id = id;
            m_has = has;
        }
 
    unsigned int masterKey() const { return m_masterKey; };
    unsigned int cfgKey() const { return m_cfgKey; };
    unsigned int cfgId() const { return m_cfgId; };
    enums::Lsf::RsdState state() const { return m_state; };
    enums::Lsf::LeakedPrescaler prescaler() const { return m_prescaler; };
    unsigned int version() const { return m_version; };
    bool has() const { return m_has; }


      GammaHandlerRsdV0(const GammaHandlerRsdV0 &other)  {
          m_type = other.m_type;
          m_masterKey = other.m_masterKey;
          m_cfgKey = other.m_cfgKey;
          m_cfgId = other.m_cfgId;
          m_state = other.m_state;
          m_prescaler = other.m_prescaler;
          m_version = other.m_version;
          m_id = other.m_id;
          m_has = other.m_has;
          m_status = other.m_status;
          m_stage = other.m_stage;
          m_energyValid = other.m_energyValid;
          m_energyInLeus = other.m_energyInLeus;
      } 

  private:
    enums::Lsf::HandlerType       m_type;
    unsigned int      m_masterKey;  /// FMX key of master CDM for handler (fixed for a run)
    unsigned int      m_cfgKey;     /// FMX key of current CDM for handler (may vary by mode)
    unsigned int      m_cfgId;      /// unique identifier of handler config (may vary by mode)
    enums::Lsf::RsdState   m_state;      /// overall filter result for the event
    enums::Lsf::LeakedPrescaler   m_prescaler;  /// Unsupported in RSD V0, otherwise indicates how event was leaked
    unsigned int     m_version;    /// Encoding version of handler-specific RSD
    enums::Lsf::HandlerId         m_id;         /// Handler identifier
    bool              m_has;        /// Indicates whether handler generated summary data
    unsigned int m_status;
    unsigned int m_stage;
    unsigned int m_energyValid;
    signed int m_energyInLeus;
  };

  class HipHandlerRsdV0  {
  public:
    HipHandlerRsdV0() { };
   ~HipHandlerRsdV0() { };

    //HipHandlerRsdV0* clone() const {
    //  std::cout<<"Hip clone" << std::endl;
    //  return new HipHandlerRsdV0(*this);
   // }
    /// This is a poor-man's dynamic cast
    //const HipHandlerRsdV0* castToHipRsdV0() const { return this; };

    void setStatus(unsigned int status) { m_status = status;}
 
    unsigned int status() const { return m_status; }
    enums::Lsf::HandlerId id() const { return m_id; }    
 
    void set(unsigned int masterKey, unsigned int cfgKey, unsigned int cfgId, 
            enums::Lsf::RsdState state, enums::Lsf::LeakedPrescaler prescaler, 
            unsigned int version, enums::Lsf::HandlerId id, bool has) {
            m_masterKey = masterKey;
            m_cfgKey = cfgKey;
            m_cfgId = cfgId;
            m_state = state;
            m_prescaler = prescaler;
            m_version = version;
            m_id = id;
            m_has = has;
        }
 
    unsigned int masterKey() const { return m_masterKey; };
    unsigned int cfgKey() const { return m_cfgKey; };
    unsigned int cfgId() const { return m_cfgId; };
    enums::Lsf::RsdState state() const { return m_state; };
    enums::Lsf::LeakedPrescaler prescaler() const { return m_prescaler; };
    unsigned int version() const { return m_version; };
    bool has() const { return m_has; }


      HipHandlerRsdV0(const HipHandlerRsdV0 &other)  {
          m_type = other.m_type;
          m_masterKey = other.m_masterKey;
          m_cfgKey = other.m_cfgKey;
          m_cfgId = other.m_cfgId;
          m_state = other.m_state;
          m_prescaler = other.m_prescaler;
          m_version = other.m_version;
          m_id = other.m_id;
          m_has = other.m_has;
          m_status = other.m_status;
      } 

  private:
    enums::Lsf::HandlerType       m_type;
    unsigned int      m_masterKey;  /// FMX key of master CDM for handler (fixed for a run)
    unsigned int      m_cfgKey;     /// FMX key of current CDM for handler (may vary by mode)
    unsigned int      m_cfgId;      /// unique identifier of handler config (may vary by mode)
    enums::Lsf::RsdState   m_state;      /// overall filter result for the event
    enums::Lsf::LeakedPrescaler   m_prescaler;  /// Unsupported in RSD V0, otherwise indicates how event was leaked
    unsigned int     m_version;    /// Encoding version of handler-specific RSD
    enums::Lsf::HandlerId         m_id;         /// Handler identifier
    bool              m_has;        /// Indicates whether handler generated summary data
    unsigned int m_status;
  };

  class MipHandlerRsdV0  {
  public:
      MipHandlerRsdV0()  { };
      ~MipHandlerRsdV0() { };

    //MipHandlerRsdV0* clone() const {
    // std::cout<<"Mip clone" << std::endl;
    //  return new MipHandlerRsdV0(*this);
   // }
    /// This is a poor-man's dynamic cast
    //const MipHandlerRsdV0* castToMipRsdV0() const { return this; };
    void setStatus(unsigned int status) { m_status = status; }

    unsigned int status() const { return m_status; }
    enums::Lsf::HandlerId id() const { return m_id; }    

    void set(unsigned int masterKey, unsigned int cfgKey, unsigned int cfgId, 
            enums::Lsf::RsdState state, enums::Lsf::LeakedPrescaler prescaler, 
            unsigned int version, enums::Lsf::HandlerId id, bool has) {
            m_masterKey = masterKey;
            m_cfgKey = cfgKey;
            m_cfgId = cfgId;
            m_state = state;
            m_prescaler = prescaler;
            m_version = version;
            m_id = id;
            m_has = has;
        }
 
    unsigned int masterKey() const { return m_masterKey; };
    unsigned int cfgKey() const { return m_cfgKey; };
    unsigned int cfgId() const { return m_cfgId; };
    enums::Lsf::RsdState state() const { return m_state; };
    enums::Lsf::LeakedPrescaler prescaler() const { return m_prescaler; };
    unsigned int version() const { return m_version; };
    bool has() const { return m_has; }


      MipHandlerRsdV0(const MipHandlerRsdV0 &other)  {
          m_type = other.m_type;
          m_masterKey = other.m_masterKey;
          m_cfgKey = other.m_cfgKey;
          m_cfgId = other.m_cfgId;
          m_state = other.m_state;
          m_prescaler = other.m_prescaler;
          m_version = other.m_version;
          m_id = other.m_id;
          m_has = other.m_has;
          m_status = other.m_status;
      } 

  private:
    enums::Lsf::HandlerType       m_type;
    unsigned int      m_masterKey;  /// FMX key of master CDM for handler (fixed for a run)
    unsigned int      m_cfgKey;     /// FMX key of current CDM for handler (may vary by mode)
    unsigned int      m_cfgId;      /// unique identifier of handler config (may vary by mode)
    enums::Lsf::RsdState   m_state;      /// overall filter result for the event
    enums::Lsf::LeakedPrescaler   m_prescaler;  /// Unsupported in RSD V0, otherwise indicates how event was leaked
    unsigned int     m_version;    /// Encoding version of handler-specific RSD
    enums::Lsf::HandlerId         m_id;         /// Handler identifier
    bool              m_has;        /// Indicates whether handler generated summary data
    unsigned int m_status;
  };


  class PassthruHandlerRsdV0  {
  public:
    PassthruHandlerRsdV0() {};
    ~PassthruHandlerRsdV0() {};

     //PassthruHandlerRsdV0* clone() const {
     // std::cout << "pass clone " << std::endl;
     // return new PassthruHandlerRsdV0(*this);
   // }
    /// This is a poor-man's dynamic cast
    //const PassthruHandlerRsdV0* castToPassthruRsdV0() const { return this; };

    void setStatus(unsigned int status) { m_status = status; }
 
    unsigned int status() const { return m_status; }
    enums::Lsf::HandlerId id() const { return m_id; }    

    void set(unsigned int masterKey, unsigned int cfgKey, unsigned int cfgId, 
            enums::Lsf::RsdState state, enums::Lsf::LeakedPrescaler prescaler, 
            unsigned int version, enums::Lsf::HandlerId id, bool has) {
            m_masterKey = masterKey;
            m_cfgKey = cfgKey;
            m_cfgId = cfgId;
            m_state = state;
            m_prescaler = prescaler;
            m_version = version;
            m_id = id;
            m_has = has;
        }
 
    unsigned int masterKey() const { return m_masterKey; };
    unsigned int cfgKey() const { return m_cfgKey; };
    unsigned int cfgId() const { return m_cfgId; };
    enums::Lsf::RsdState state() const { return m_state; };
    enums::Lsf::LeakedPrescaler prescaler() const { return m_prescaler; };
    unsigned int version() const { return m_version; };
    bool has() const { return m_has; }


      PassthruHandlerRsdV0(const PassthruHandlerRsdV0 &other)  {
          m_type = other.m_type;
          m_masterKey = other.m_masterKey;
          m_cfgKey = other.m_cfgKey;
          m_cfgId = other.m_cfgId;
          m_state = other.m_state;
          m_prescaler = other.m_prescaler;
          m_version = other.m_version;
          m_id = other.m_id;
          m_has = other.m_has;
          m_status = other.m_status;
      } 

  private:
    enums::Lsf::HandlerType       m_type;
    unsigned int      m_masterKey;  /// FMX key of master CDM for handler (fixed for a run)
    unsigned int      m_cfgKey;     /// FMX key of current CDM for handler (may vary by mode)
    unsigned int      m_cfgId;      /// unique identifier of handler config (may vary by mode)
    enums::Lsf::RsdState   m_state;      /// overall filter result for the event
    enums::Lsf::LeakedPrescaler   m_prescaler;  /// Unsupported in RSD V0, otherwise indicates how event was leaked
    unsigned int     m_version;    /// Encoding version of handler-specific RSD
    enums::Lsf::HandlerId         m_id;         /// Handler identifier
    bool              m_has;        /// Indicates whether handler generated summary data
    unsigned int m_status;
  };


}

#endif // LSFDATA_LPAHANDLER_HH

