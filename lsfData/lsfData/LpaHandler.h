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


    enums::Lsf::HandlerId id() const { return m_id; }    


    void print(const std::string &str="") const { }

    friend std::ostream& operator<< ( std::ostream& s, const LpaHandler& o) { return s; };
    

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

    void setStatus(unsigned int status) { m_status = status; }
    unsigned int status() const { return m_status; }

      DgnHandlerRsdV0(const DgnHandlerRsdV0 &other)  {
          m_status = other.m_status;
      } 

  private:
    unsigned int m_status;


  };
    
class DgnHandler {
public:
    DgnHandler():m_dgn(0) { };

    DgnHandler(const DgnHandler &other)  {
          m_handler = other.m_handler;
          if (other.m_dgn)
              m_dgn = new DgnHandlerRsdV0(*(other.m_dgn));
           else
              m_dgn = 0;
      } 

   // const char*                     typeName() const;
    ~DgnHandler() { 
       if (m_dgn) {
           delete m_dgn;
           m_dgn = 0;
        }
     };

    void set(unsigned int masterKey, unsigned int cfgKey, unsigned int cfgId, 
            enums::Lsf::RsdState state, enums::Lsf::LeakedPrescaler prescaler, 
            unsigned int version, enums::Lsf::HandlerId id, bool has) {
        m_handler.set(masterKey, cfgKey, cfgId, state, prescaler,
             version, id, has);
    }

    void setLpaHandler(const LpaHandler &lpaHandler) {
        m_handler = lpaHandler;
    }

    //void setRsd(const DgnHandlerRsdV0* dgn) { m_dgn = dgn; }

    void setStatus(unsigned int status) {
        if (!m_dgn) m_dgn = new DgnHandlerRsdV0;
        m_dgn->setStatus(status);
    }

    const LpaHandler& lpaHandler() const { return m_handler; }

    const DgnHandlerRsdV0* rsd() const { return m_dgn; }
    unsigned int masterKey() const { return m_handler.masterKey(); };
    unsigned int cfgKey() const { return m_handler.cfgKey(); };
    unsigned int cfgId() const { return m_handler.cfgId(); };
    enums::Lsf::RsdState state() const { return m_handler.state(); };
    enums::Lsf::LeakedPrescaler prescaler() const { return m_handler.prescaler(); };
    unsigned int version() const { return m_handler.version(); };
    bool has() const { return m_handler.has(); }
    enums::Lsf::HandlerId id() const { return m_handler.id(); }    


private:
  LpaHandler m_handler;
  DgnHandlerRsdV0 *m_dgn;


};


  class GammaHandlerRsdV0  {
  public:
      GammaHandlerRsdV0()  { };
      ~GammaHandlerRsdV0() { };

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
 
      GammaHandlerRsdV0(const GammaHandlerRsdV0 &other)  {
          m_status = other.m_status;
          m_stage = other.m_stage;
          m_energyValid = other.m_energyValid;
          m_energyInLeus = other.m_energyInLeus;
      } 

  private:
    unsigned int m_status;
    unsigned int m_stage;
    unsigned int m_energyValid;
    signed int m_energyInLeus;
  };

class GammaHandler {
public:
    GammaHandler():m_gamma(0) { };
    GammaHandler(const GammaHandler &other)  {
          m_handler = other.m_handler;
          if (other.m_gamma)
              m_gamma = new GammaHandlerRsdV0(*(other.m_gamma));
          else 
              m_gamma = 0;
      } 

    ~GammaHandler() { 
       if (m_gamma) {
           delete m_gamma;
           m_gamma = 0;
        }
     };

    void set(unsigned int masterKey, unsigned int cfgKey, unsigned int cfgId, 
            enums::Lsf::RsdState state, enums::Lsf::LeakedPrescaler prescaler, 
            unsigned int version, enums::Lsf::HandlerId id, bool has) {
        m_handler.set(masterKey, cfgKey, cfgId, state, prescaler,
              version, id, has);
    }
    void setLpaHandler(const LpaHandler &lpaHandler) {
        m_handler = lpaHandler;
    }

    void setStatus(unsigned int status, unsigned int stage,
                   unsigned int energyValid, int energyInLeus) {
            if (!m_gamma) m_gamma = new GammaHandlerRsdV0;
            m_gamma->setStatus(status, stage, energyValid, energyInLeus);
        }
    //void setRsd(const GammaHandlerRsdV0* gam) { m_gamma = gam; }

    const LpaHandler& lpaHandler() const { return m_handler; }

    const GammaHandlerRsdV0* rsd() const { return m_gamma; }
    unsigned int masterKey() const { return m_handler.masterKey(); };
    unsigned int cfgKey() const { return m_handler.cfgKey(); };
    unsigned int cfgId() const { return m_handler.cfgId(); };
    enums::Lsf::RsdState state() const { return m_handler.state(); };
    enums::Lsf::LeakedPrescaler prescaler() const { return m_handler.prescaler(); };
    unsigned int version() const { return m_handler.version(); };
    bool has() const { return m_handler.has(); }
    enums::Lsf::HandlerId id() const { return m_handler.id(); }    

private:
  LpaHandler m_handler;
  GammaHandlerRsdV0 *m_gamma;


};


  class HipHandlerRsdV0  {
  public:
    HipHandlerRsdV0() { };
   ~HipHandlerRsdV0() { };

    void setStatus(unsigned int status) { m_status = status;}
 
    unsigned int status() const { return m_status; }
 
      HipHandlerRsdV0(const HipHandlerRsdV0 &other)  {
          m_status = other.m_status;
      } 

  private:
    unsigned int m_status;
  };

class HipHandler {
public:
    HipHandler():m_hip(0) { };

    HipHandler(const HipHandler &other)  {
          m_handler = other.m_handler;
          if (other.m_hip)
              m_hip = new HipHandlerRsdV0(*(other.m_hip));
          else
              m_hip = 0;
      } 
    ~HipHandler() { 
       if (m_hip) {
           delete m_hip;
           m_hip = 0;
        }
     };

    void set(unsigned int masterKey, unsigned int cfgKey, unsigned int cfgId, 
            enums::Lsf::RsdState state, enums::Lsf::LeakedPrescaler prescaler, 
            unsigned int version, enums::Lsf::HandlerId id, bool has) {
        m_handler.set(masterKey, cfgKey, cfgId, state, prescaler,
             version, id, has);
    }
    void setLpaHandler(const LpaHandler &lpaHandler) {
        m_handler = lpaHandler;
    }

    void setStatus(unsigned int status) {
        if (!m_hip) m_hip = new HipHandlerRsdV0;
        m_hip->setStatus(status);
    }
    //void setRsd(const HipHandlerRsdV0* hip) { m_hip = hip; }

    const LpaHandler& lpaHandler() const { return m_handler; }

    const HipHandlerRsdV0* rsd() const { return m_hip; }
    unsigned int masterKey() const { return m_handler.masterKey(); };
    unsigned int cfgKey() const { return m_handler.cfgKey(); };
    unsigned int cfgId() const { return m_handler.cfgId(); };
    enums::Lsf::RsdState state() const { return m_handler.state(); };
    enums::Lsf::LeakedPrescaler prescaler() const { return m_handler.prescaler(); };
    unsigned int version() const { return m_handler.version(); };
    bool has() const { return m_handler.has(); }
    enums::Lsf::HandlerId id() const { return m_handler.id(); }    

private:
  LpaHandler m_handler;
  HipHandlerRsdV0 *m_hip;


};


  class MipHandlerRsdV0  {
  public:
      MipHandlerRsdV0()  { };
      ~MipHandlerRsdV0() { };

    void setStatus(unsigned int status) { m_status = status; }

    unsigned int status() const { return m_status; }

      MipHandlerRsdV0(const MipHandlerRsdV0 &other)  {
          m_status = other.m_status;
      } 

  private:
    unsigned int m_status;
  };


class MipHandler {
public:
    MipHandler():m_mip(0) { };

    MipHandler(const MipHandler &other)  {
          m_handler = other.m_handler;
          if (other.m_mip)
              m_mip = new MipHandlerRsdV0(*(other.m_mip));
          else
              m_mip = 0;
      } 
    ~MipHandler() { 
       if (m_mip) {
           delete m_mip;
           m_mip = 0;
        }
     };

    void set(unsigned int masterKey, unsigned int cfgKey, unsigned int cfgId, 
            enums::Lsf::RsdState state, enums::Lsf::LeakedPrescaler prescaler, 
            unsigned int version, enums::Lsf::HandlerId id, bool has) {
        m_handler.set(masterKey, cfgKey, cfgId, state, prescaler,
             version, id, has);
    }
    void setLpaHandler(const LpaHandler &lpaHandler) {
        m_handler = lpaHandler;
    }

    //void setRsd(const MipHandlerRsdV0* mip) { m_mip = mip; }

    const LpaHandler& lpaHandler() const { return m_handler; }

    void setStatus(unsigned int status) {
        if (!m_mip) m_mip = new MipHandlerRsdV0;
        m_mip->setStatus(status);
    }
    const MipHandlerRsdV0* rsd() const { return m_mip; }
    unsigned int masterKey() const { return m_handler.masterKey(); };
    unsigned int cfgKey() const { return m_handler.cfgKey(); };
    unsigned int cfgId() const { return m_handler.cfgId(); };
    enums::Lsf::RsdState state() const { return m_handler.state(); };
    enums::Lsf::LeakedPrescaler prescaler() const { return m_handler.prescaler(); };
    unsigned int version() const { return m_handler.version(); };
    bool has() const { return m_handler.has(); }
    enums::Lsf::HandlerId id() const { return m_handler.id(); }    

private:
  LpaHandler m_handler;
  MipHandlerRsdV0 *m_mip;



};


  class PassthruHandlerRsdV0  {
  public:
    PassthruHandlerRsdV0() {};
    ~PassthruHandlerRsdV0() {};

    void setStatus(unsigned int status) { m_status = status; }
 
    unsigned int status() const { return m_status; }

      PassthruHandlerRsdV0(const PassthruHandlerRsdV0 &other)  {
          m_status = other.m_status;
      } 

  private:
    unsigned int m_status;
  };

class PassthruHandler {
public:
    PassthruHandler():m_pass(0) { };

    PassthruHandler(const PassthruHandler &other)  {
          m_handler = other.m_handler;
          if (other.m_pass)
              m_pass = new PassthruHandlerRsdV0(*(other.m_pass));
          else
              m_pass = 0;
      } 
    ~PassthruHandler() { 
       if (m_pass) {
           delete m_pass;
           m_pass = 0;
        }
     };

    void set(unsigned int masterKey, unsigned int cfgKey, unsigned int cfgId, 
            enums::Lsf::RsdState state, enums::Lsf::LeakedPrescaler prescaler, 
            unsigned int version, enums::Lsf::HandlerId id, bool has) {
        m_handler.set(masterKey, cfgKey, cfgId, state, prescaler,
             version, id, has);
    }
    void setLpaHandler(const LpaHandler &lpaHandler) {
        m_handler = lpaHandler;
    }

    void setStatus(unsigned int status) {
        if (!m_pass) m_pass = new PassthruHandlerRsdV0;
        m_pass->setStatus(status);
    }
    //void setRsd(const PassthruHandlerRsdV0* pass) { m_pass = pass; }

    const LpaHandler& lpaHandler() const { return m_handler; }

    const PassthruHandlerRsdV0* rsd() const { return m_pass; }
    unsigned int masterKey() const { return m_handler.masterKey(); };
    unsigned int cfgKey() const { return m_handler.cfgKey(); };
    unsigned int cfgId() const { return m_handler.cfgId(); };
    enums::Lsf::RsdState state() const { return m_handler.state(); };
    enums::Lsf::LeakedPrescaler prescaler() const { return m_handler.prescaler(); };
    unsigned int version() const { return m_handler.version(); };
    bool has() const { return m_handler.has(); }
    enums::Lsf::HandlerId id() const { return m_handler.id(); }    

private:
  LpaHandler m_handler;
  PassthruHandlerRsdV0 *m_pass;



};


}

#endif // LSFDATA_LPAHANDLER_HH

