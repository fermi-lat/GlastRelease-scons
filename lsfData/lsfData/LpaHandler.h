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
 class DgnRsdV0;
 class GammaRsdV0;
 class GammaRsdV1;
 class GammaRsdV2;
 class GammaRsdV3;
 class HipRsdV0;
 class MipRsdV0;
 class PassthruRsdV0;


  /** class representing the result info from each active handler */
  class LpaHandler { 
  public:

    /** ctor and dump routine */
      LpaHandler():m_type(enums::Lsf::Unknown), m_masterKey(0xFFFFFFFF), m_cfgKey(0xFFFFFFFF), 
          m_cfgId(0xFFFFFFFF), m_state(enums::Lsf::INVALID), m_prescaler(enums::Lsf::UNSUPPORTED), m_version(0), 
          m_id(enums::Lsf::HandlerIdCnt), m_has(false), 
          m_prescaleFactor(LSF_INVALID_UINT) {};

      virtual ~LpaHandler() { };


    enums::Lsf::HandlerId id() const { return m_id; }    


    void print(const std::string &str="") const { printf("%s\n", str);}

    friend std::ostream& operator<< ( std::ostream& s, const LpaHandler& /*o*/) { return s; };
    

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

    void setPrescaleFactor(unsigned int prescaleFactor) {
        m_prescaleFactor = prescaleFactor; 
    }
 
    unsigned int masterKey() const { return m_masterKey; };
    unsigned int cfgKey() const { return m_cfgKey; };
    unsigned int cfgId() const { return m_cfgId; };
    enums::Lsf::RsdState state() const { return m_state; };
    enums::Lsf::LeakedPrescaler prescaler() const { return m_prescaler; };
    unsigned int version() const { return m_version; };
    bool has() const { return m_has; }
    unsigned int prescaleFactor() const { return m_prescaleFactor; };

    int prescaleIndex() const {
      switch (m_prescaler){
      case enums::Lsf::UNSUPPORTED:
	// We really don't know what happened, maybe the state will tell us something
	break;
      case enums::Lsf::INPUT:
	// filter never saw the event, set a value larger than any specific prescaler
	return 33;
      case enums::Lsf::OUTPUT:
	// output prescale expired, this means the line prescalers are irrelevent
	// set a value larger that the first line
	return 32;      
      default:
	// one of the line prescales is asserted, return that 
	return m_prescaler;
      }
      switch (m_state) {
      case enums::Lsf::INVALID:  
	// We really don't know what happened, set a large value so this 
	// doesn't get into any of the other samples
	return 100;
      case enums::Lsf::IGNORED:
	// Filter never saw the event, set a value larger than any specific prescaler
	// should have been caught in the first switch above
	return 33;
      case enums::Lsf::PASSED:
	// If the event made it here it means that none of the bits were asserted
	// return a value smaller that the lowest bit
	return -1;
      case enums::Lsf::SUPPRESSED:
      case enums::Lsf::VETOED:
      case enums::Lsf::LEAKED:  
	// In all of these cases the event should have been caught by the first switch
	// but return the relevent prescaler, just in case
	break;
      }
      return m_prescaler;
    }



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
          m_prescaleFactor = other.m_prescaleFactor;
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
    unsigned int     m_prescaleFactor;


  };

  /** classes representing V0 of the handler-specific RSD's */
  class DgnRsdV0  {
  public:
      DgnRsdV0() {};
      ~DgnRsdV0() {};

    void setStatus(unsigned int status) { m_status = status; }
    unsigned int status() const { return m_status; }

    unsigned int stage() const {
        return m_status & DgnRsdV0::Stage_m;
    }

    // status-word value accessor
    unsigned int gem_classes() const {
      return (m_status & DgnRsdV0::Gem_Classes_m) >> DgnRsdV0::Gem_Classes_v;
  }

      DgnRsdV0(const DgnRsdV0 &other)  {
          m_status = other.m_status;
      } 

  private:
    unsigned int m_status;

   enum Status_s {
      Gem_Classes_s = 16,
    };
    enum Status_v {
      Stage_Gem_v   = 0,
      Gem_Classes_v = 15,
    };
    enum Status_m {
      Stage_Gem_m   = 1 << Stage_Gem_v,
      Gem_Classes_m = ((1 << Gem_Classes_s) - 1) << Gem_Classes_v,
      Stage_m       = Stage_Gem_m,
    };


  };
    
class DgnHandler {
public:
    DgnHandler():m_dgn(0) { };

    DgnHandler(const DgnHandler &other)  {
          m_handler = other.m_handler;
          if (other.m_dgn)
              m_dgn = new DgnRsdV0(*(other.m_dgn));
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
            unsigned int version, enums::Lsf::HandlerId id, bool has,
            unsigned int prescaleFactor=LSF_INVALID_UINT) {
        m_handler.set(masterKey, cfgKey, cfgId, state, prescaler,
             version, id, has);
        m_handler.setPrescaleFactor(prescaleFactor);
    }

    void setLpaHandler(const LpaHandler &lpaHandler) {
        m_handler = lpaHandler;
    }

    void setPrescaleFactor(unsigned int factor) { 
        m_handler.setPrescaleFactor(factor); 
    }


    void setStatus(unsigned int status) {
        if (!m_dgn) m_dgn = new DgnRsdV0;
        m_dgn->setStatus(status);
    }

    const LpaHandler& lpaHandler() const { return m_handler; }

    const DgnRsdV0* rsd() const { return m_dgn; }
    unsigned int masterKey() const { return m_handler.masterKey(); };
    unsigned int cfgKey() const { return m_handler.cfgKey(); };
    unsigned int cfgId() const { return m_handler.cfgId(); };
    enums::Lsf::RsdState state() const { return m_handler.state(); };
    enums::Lsf::LeakedPrescaler prescaler() const { return m_handler.prescaler(); };
    unsigned int version() const { return m_handler.version(); };
    bool has() const { return m_handler.has(); }
    enums::Lsf::HandlerId id() const { return m_handler.id(); }    
    unsigned int prescaleFactor() const { return m_handler.prescaleFactor(); }


private:
  LpaHandler m_handler;
  DgnRsdV0 *m_dgn;


};

  class GammaRsd  {
  public:
      GammaRsd()  { };
      virtual ~GammaRsd() { };

    void setStatus(unsigned int status, unsigned int stage, unsigned int energyValid, int energyInLeus) {
            m_status = status;
            m_stage = stage;
            m_energyValid = energyValid;
            m_energyInLeus = energyInLeus;
        }
 
     GammaRsd(const GammaRsd &other)  {
          m_status = other.m_status;
          m_stage = other.m_stage;
          m_energyValid = other.m_energyValid;
          m_energyInLeus = other.m_energyInLeus;
      } 

    virtual unsigned int status() const = 0;
    virtual unsigned int stage() const = 0;
    virtual unsigned int energyValid() const = 0; 
    virtual int energyInLeus() const = 0;

  protected:
    unsigned int m_status;
    unsigned int m_stage;
    unsigned int m_energyValid;
    signed int m_energyInLeus;
  };


  class GammaRsdV0 : public GammaRsd  {
  public:
      GammaRsdV0():GammaRsd()  { };
      virtual ~GammaRsdV0() { };

    unsigned int status() const { return m_status; }
    unsigned int stage() const { return m_stage; }
    unsigned int energyValid() const { return m_energyValid; }
    int energyInLeus() const { return m_energyInLeus; }

      GammaRsdV0(const GammaRsdV0 &other):GammaRsd(other)  {
          m_status = other.m_status;
          m_stage = other.m_stage;
          m_energyValid = other.m_energyValid;
          m_energyInLeus = other.m_energyInLeus;
      } 

      GammaRsdV0(const GammaRsd &other) : GammaRsd(other) { }
          

  };

  class GammaRsdV1 : public GammaRsd {
  public:
      GammaRsdV1():GammaRsd()  { };
      virtual ~GammaRsdV1() { };

    unsigned int status() const { return m_status; }
    unsigned int stage() const { return m_stage; }
    unsigned int energyValid() const { return m_energyValid; }
    int energyInLeus() const { return m_energyInLeus; }
 
      GammaRsdV1(const GammaRsdV1 &other):GammaRsd(other)  {
          m_status = other.m_status;
          m_stage = other.m_stage;
          m_energyValid = other.m_energyValid;
          m_energyInLeus = other.m_energyInLeus;
      } 

      GammaRsdV1(const GammaRsd &other) : GammaRsd(other) { }

  };

  class GammaRsdV2 : public GammaRsd {
  public:
      GammaRsdV2():GammaRsd()  { };
      virtual ~GammaRsdV2() { };

    unsigned int status() const { return m_status; }
    unsigned int stage() const { return m_stage; }
    unsigned int energyValid() const { return m_energyValid; }
    int energyInLeus() const { return m_energyInLeus; }
 
      GammaRsdV2(const GammaRsdV2 &other):GammaRsd(other)  {
          m_status = other.m_status;
          m_stage = other.m_stage;
          m_energyValid = other.m_energyValid;
          m_energyInLeus = other.m_energyInLeus;
      } 

      GammaRsdV2(const GammaRsd &other) : GammaRsd(other) { }
  };

  class GammaRsdV3 : public GammaRsd {
  public:
      GammaRsdV3():GammaRsd()  { };
      virtual ~GammaRsdV3() { };

    unsigned int status() const { return m_status; }
    unsigned int stage() const { return m_stage; }
    unsigned int energyValid() const { return m_energyValid; }
    int energyInLeus() const { return m_energyInLeus; }
 
      GammaRsdV3(const GammaRsdV3 &other):GammaRsd(other)  {
          m_status = other.m_status;
          m_stage = other.m_stage;
          m_energyValid = other.m_energyValid;
          m_energyInLeus = other.m_energyInLeus;
      } 

      GammaRsdV3(const GammaRsd &other) : GammaRsd(other) { }
  };

class GammaHandler {
public:
    GammaHandler():m_gamma(0) { };
    GammaHandler(const GammaHandler &other)  {
          m_handler = other.m_handler;
          if (other.m_gamma)
            switch(m_handler.version()) {
              case 0:
                m_gamma = new GammaRsdV0(*(other.m_gamma));
                break;
              case 1:
                m_gamma = new GammaRsdV1(*(other.m_gamma));
                break;
              case 2:
                m_gamma = new GammaRsdV2(*(other.m_gamma));
                break;
              case 3:
                m_gamma = new GammaRsdV3(*(other.m_gamma));
                break;
              default:
                std:: cout << "No valid version for GammaRsd found, set to NULL"
                           << std::endl;
                m_gamma = 0;
            }
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
            unsigned int version, enums::Lsf::HandlerId id, bool has,
            unsigned int prescaleFactor=LSF_INVALID_UINT) {
        m_handler.set(masterKey, cfgKey, cfgId, state, prescaler,
              version, id, has);
        m_handler.setPrescaleFactor(prescaleFactor); 
    }
    void setLpaHandler(const LpaHandler &lpaHandler) {
        m_handler = lpaHandler;
    }

    void setPrescaleFactor(unsigned int factor) { 
        m_handler.setPrescaleFactor(factor); 
    }
    void setStatus(unsigned int status, unsigned int stage,
                   unsigned int energyValid, int energyInLeus) {
            if (!m_gamma)
              switch(m_handler.version()) {
              case 0:
                m_gamma = new GammaRsdV0;
                break;
              case 1:
                m_gamma = new GammaRsdV1;
                break;
              case 2:
                m_gamma = new GammaRsdV2;
                break;
              case 3:
                m_gamma = new GammaRsdV3;
                break;
              default:
                m_gamma = 0;
                std::cout << "Gamma version invalid, not setting GammaRsd"
                          << std::endl;
                return;
              }
            m_gamma->setStatus(status, stage, energyValid, energyInLeus);
        }

    const LpaHandler& lpaHandler() const { return m_handler; }

    const GammaRsd* rsd() const { return m_gamma; }
    unsigned int masterKey() const { return m_handler.masterKey(); };
    unsigned int cfgKey() const { return m_handler.cfgKey(); };
    unsigned int cfgId() const { return m_handler.cfgId(); };
    enums::Lsf::RsdState state() const { return m_handler.state(); };
    enums::Lsf::LeakedPrescaler prescaler() const { return m_handler.prescaler(); };
    unsigned int version() const { return m_handler.version(); };
    bool has() const { return m_handler.has(); }
    enums::Lsf::HandlerId id() const { return m_handler.id(); }    
    unsigned int prescaleFactor() const { return m_handler.prescaleFactor(); }

private:
  LpaHandler m_handler;
  GammaRsd *m_gamma;


};


  class HipRsdV0  {
  public:
    HipRsdV0() { };
   ~HipRsdV0() { };

    void setStatus(unsigned int status) { m_status = status;}
 
    unsigned int status() const { return m_status; }
    unsigned int stage() const {
        return m_status & HipRsdV0::Stage_m;
    }
 
      HipRsdV0(const HipRsdV0 &other)  {
          m_status = other.m_status;
      } 

  private:
    unsigned int m_status;

    enum Status_v {
      Stage_Gem_v      = 0,
      Stage_Dir_v      = 1,
      Stage_Cal_v      = 2,
      Stage_Cal_Echk_v = 3,
      Stage_Cal_Lchk_v = 4,
    };
    enum Status_m {
      Stage_Gem_m      = 1 << Stage_Gem_v,
      Stage_Dir_m      = 1 << Stage_Dir_v,
      Stage_Cal_m      = 1 << Stage_Cal_v,
      Stage_Cal_Echk_m = 1 << Stage_Cal_Echk_v,
      Stage_Cal_Lchk_m = 1 << Stage_Cal_Lchk_v,
      Stage_m          = Stage_Gem_m
                       | Stage_Dir_m
                       | Stage_Cal_m
                       | Stage_Cal_Echk_m
                       | Stage_Cal_Lchk_m,
    };

  };

class HipHandler {
public:
    HipHandler():m_hip(0) { };

    HipHandler(const HipHandler &other)  {
          m_handler = other.m_handler;
          if (other.m_hip)
              m_hip = new HipRsdV0(*(other.m_hip));
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
            unsigned int version, enums::Lsf::HandlerId id, bool has,
            unsigned int prescaleFactor = LSF_INVALID_UINT) {
        m_handler.set(masterKey, cfgKey, cfgId, state, prescaler,
             version, id, has);
        m_handler.setPrescaleFactor(prescaleFactor); 
    }
    void setLpaHandler(const LpaHandler &lpaHandler) {
        m_handler = lpaHandler;
    }
    void setPrescaleFactor(unsigned int factor) { 
        m_handler.setPrescaleFactor(factor); 
    }

    void setStatus(unsigned int status) {
        if (!m_hip) m_hip = new HipRsdV0;
        m_hip->setStatus(status);
    }
    //void setRsd(const HipRsdV0* hip) { m_hip = hip; }

    const LpaHandler& lpaHandler() const { return m_handler; }

    const HipRsdV0* rsd() const { return m_hip; }
    unsigned int masterKey() const { return m_handler.masterKey(); };
    unsigned int cfgKey() const { return m_handler.cfgKey(); };
    unsigned int cfgId() const { return m_handler.cfgId(); };
    enums::Lsf::RsdState state() const { return m_handler.state(); };
    enums::Lsf::LeakedPrescaler prescaler() const { return m_handler.prescaler(); };
    unsigned int version() const { return m_handler.version(); };
    bool has() const { return m_handler.has(); }
    enums::Lsf::HandlerId id() const { return m_handler.id(); }    
    unsigned int prescaleFactor() const { return m_handler.prescaleFactor(); }

private:
  LpaHandler m_handler;
  HipRsdV0 *m_hip;


};


  class MipRsdV0  {
  public:
      MipRsdV0()  { };
      ~MipRsdV0() { };

    void setStatus(unsigned int status) { m_status = status; }

    unsigned int status() const { return m_status; }
    unsigned int stage() const {
      return m_status & MipRsdV0::Stage_m;
    }

      MipRsdV0(const MipRsdV0 &other)  {
          m_status = other.m_status;
      } 

  private:
    unsigned int m_status;

    enum Status_v {
      Stage_Gem_v   = 0,
      Stage_Acd_v   = 1,
      Stage_Dir_v   = 2,
      Stage_Cal_v   = 3,
      Stage_XCal_v  = 4,
    };
    enum Status_m {
      Stage_Gem_m   = 1 << Stage_Gem_v,
      Stage_Acd_m   = 1 << Stage_Acd_v,
      Stage_Dir_m   = 1 << Stage_Dir_v,
      Stage_Cal_m   = 1 << Stage_Cal_v,
      Stage_XCal_m  = 1 << Stage_XCal_v,
      Stage_m       = Stage_Gem_m
                    | Stage_Acd_m
                    | Stage_Dir_m
                    | Stage_Cal_m
                    | Stage_XCal_m,
    };
  };


class MipHandler {
public:
    MipHandler():m_mip(0) { };

    MipHandler(const MipHandler &other)  {
          m_handler = other.m_handler;
          if (other.m_mip)
              m_mip = new MipRsdV0(*(other.m_mip));
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
            unsigned int version, enums::Lsf::HandlerId id, bool has,
            unsigned int prescaleFactor=LSF_INVALID_UINT) {
        m_handler.set(masterKey, cfgKey, cfgId, state, prescaler,
             version, id, has);
        m_handler.setPrescaleFactor(prescaleFactor); 
    }
    void setLpaHandler(const LpaHandler &lpaHandler) {
        m_handler = lpaHandler;
    }
    void setPrescaleFactor(unsigned int factor) { 
        m_handler.setPrescaleFactor(factor); 
    }

    //void setRsd(const MipRsdV0* mip) { m_mip = mip; }

    const LpaHandler& lpaHandler() const { return m_handler; }

    void setStatus(unsigned int status) {
        if (!m_mip) m_mip = new MipRsdV0;
        m_mip->setStatus(status);
    }
    const MipRsdV0* rsd() const { return m_mip; }
    unsigned int masterKey() const { return m_handler.masterKey(); };
    unsigned int cfgKey() const { return m_handler.cfgKey(); };
    unsigned int cfgId() const { return m_handler.cfgId(); };
    enums::Lsf::RsdState state() const { return m_handler.state(); };
    enums::Lsf::LeakedPrescaler prescaler() const { return m_handler.prescaler(); };
    unsigned int version() const { return m_handler.version(); };
    bool has() const { return m_handler.has(); }
    enums::Lsf::HandlerId id() const { return m_handler.id(); }    
    unsigned int prescaleFactor() const { return m_handler.prescaleFactor(); }

private:
  LpaHandler m_handler;
  MipRsdV0 *m_mip;

};


  class PassthruRsdV0  {
  public:
    PassthruRsdV0() {};
    ~PassthruRsdV0() {};

    void setStatus(unsigned int status) { m_status = status; }
 
    unsigned int status() const { return m_status; }
    unsigned int stage() const { return 0; };

      PassthruRsdV0(const PassthruRsdV0 &other)  {
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
              m_pass = new PassthruRsdV0(*(other.m_pass));
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
            unsigned int version, enums::Lsf::HandlerId id, bool has,
            unsigned int prescaleFactor=LSF_INVALID_UINT) {
        m_handler.set(masterKey, cfgKey, cfgId, state, prescaler,
             version, id, has);
        m_handler.setPrescaleFactor(prescaleFactor); 
    }
    void setLpaHandler(const LpaHandler &lpaHandler) {
        m_handler = lpaHandler;
    }
    void setPrescaleFactor(unsigned int factor) { 
        m_handler.setPrescaleFactor(factor); 
    }

    void setStatus(unsigned int status) {
        if (!m_pass) m_pass = new PassthruRsdV0;
        m_pass->setStatus(status);
    }
    //void setRsd(const PassthruRsdV0* pass) { m_pass = pass; }

    const LpaHandler& lpaHandler() const { return m_handler; }

    const PassthruRsdV0* rsd() const { return m_pass; }
    unsigned int masterKey() const { return m_handler.masterKey(); };
    unsigned int cfgKey() const { return m_handler.cfgKey(); };
    unsigned int cfgId() const { return m_handler.cfgId(); };
    enums::Lsf::RsdState state() const { return m_handler.state(); };
    enums::Lsf::LeakedPrescaler prescaler() const { return m_handler.prescaler(); };
    unsigned int version() const { return m_handler.version(); };
    bool has() const { return m_handler.has(); }
    enums::Lsf::HandlerId id() const { return m_handler.id(); }    
    unsigned int prescaleFactor() const { return m_handler.prescaleFactor(); }

private:
  LpaHandler m_handler;
  PassthruRsdV0 *m_pass;



};


}

#endif // LSFDATA_LPAHANDLER_HH

