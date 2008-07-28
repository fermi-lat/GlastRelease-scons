#ifndef LpaHandler_H
#define LpaHandler_H

#include "TObject.h"
#include "TObjArray.h"
#include "enums/Lsf.h"
#include "enums/EfcGammaStatus.h"
#include "enums/XfcHFCStatus.h"
#include "enums/XfcMipStatus.h"
#include "enums/XfcDFCStatus.h"



/** @class ILpaHandler
* @brief Interface/base class which describes access to the results of the various 
*        onboard filters. Currently, all filters output a 32 bit status word but then
*        can also include some number of other items. This interface class deals 
*        entirely with the status bits
*
*/
class ILpaHandler : public TObject
{
public:
    ILpaHandler() : TObject(), 
                 m_masterKey(0xFFFFFFFF), m_cfgKey(0xFFFFFFFF), 
                 m_cfgId(0xFFFFFFFF), m_version(0), m_prescaleFactor(LSF_INVALID_UINT), 
                 m_state(enums::Lsf::INVALID), 
                 m_prescaler(enums::Lsf::UNSUPPORTED), 
                 m_type(enums::Lsf::Unknown),
                 m_id(enums::Lsf::MaxHandlerIds), m_has(kFALSE) {} ;

    virtual ~ILpaHandler() {}

    /// TObject functions
    virtual void Clear(Option_t *option ="");
    virtual void Print(Option_t *option="") const;

    void initialize(UInt_t masterKey, UInt_t cfgKey, UInt_t cfgId, 
                    enums::Lsf::RsdState state, 
                    enums::Lsf::LeakedPrescaler prescaler, UInt_t version, 
                    enums::Lsf::HandlerId id, Bool_t has,
                    UInt_t prescaleFactor=LSF_INVALID_UINT);

    void setPrescaleFactor(UInt_t p) { m_prescaleFactor = p; }


    /// Lpa functions
    virtual UInt_t  getStatusWord() const { return 0; };

    virtual Bool_t passed() const { 
        if (m_state == enums::Lsf::PASSED) return true;
        else return false;
    }

    virtual Bool_t vetoed() const {
        if(m_state == enums::Lsf::VETOED) return true;
        else return false;
    }

    virtual Bool_t leaked() const {
        if (m_state == enums::Lsf::LEAKED) return true;
        else return false;
    }

    virtual Bool_t suppressed() const {
        if (m_state == enums::Lsf::SUPPRESSED) return true;
        else return false;
    }

    virtual Bool_t ignored() const {
        if (m_state == enums::Lsf::IGNORED) return true;
        return false;
   }

    UInt_t getMasterKey() const { return m_masterKey; };
    UInt_t getCfgKey() const { return m_cfgKey; };
    UInt_t getCfgId() const { return m_cfgId; };
    enums::Lsf::RsdState getState() const { return m_state; };
    enums::Lsf::LeakedPrescaler getPrescaler() const { return m_prescaler; };
    UInt_t getVersion() const { return m_version; };
    Bool_t has() const { return m_has; }
    enums::Lsf::HandlerId getId() const { return m_id; }
    UInt_t getPrescaleFactor() const { return m_prescaleFactor; }

    Int_t prescalerIndex() const;

private:
    UInt_t                      m_masterKey;  /// FMX key of master CDM for handler (fixed for a run)
    UInt_t                      m_cfgKey;     /// FMX key of current CDM for handler (may vary by mode)
    UInt_t                      m_cfgId;      /// unique identifier of handler config (may vary by mode)
    UInt_t                      m_version;    /// Encoding version of handler-specific RSD
    UInt_t                      m_prescaleFactor;
    enums::Lsf::RsdState        m_state;      /// overall filter result for the event
    enums::Lsf::LeakedPrescaler m_prescaler;  /// Unsupported in RSD V0, otherwise indicates how event was leaked
    enums::Lsf::HandlerType     m_type;
    enums::Lsf::HandlerId       m_id;         /// Handler identifier
    Bool_t                      m_has;        /// Indicates whether handler generated summary data

    ClassDef(ILpaHandler,1) // ILpaHandler interface class for all LpaHandlers

};


class LpaGammaFilter : public ILpaHandler
{
public:
    LpaGammaFilter() : m_status(0) {}
    LpaGammaFilter(UInt_t status) : m_status(status) {}
    virtual ~LpaGammaFilter() {}

    virtual void set(UInt_t status, UInt_t stage, UInt_t energyValid, Int_t energyInLeus) {
        m_status = status;
        m_stage = stage;
        m_energyValid = energyValid;
        m_energyInLeus = energyInLeus;
    }

    virtual UInt_t getStatusWord() const { return m_status; };
    virtual UInt_t getStage() const { return m_stage; };
    virtual UInt_t getEnergyValid() const { return m_energyValid; }
    virtual Int_t getEnergyInLeus() const { return m_energyInLeus; }

    virtual UInt_t getAllVetoBits() const;

    virtual void Clear(Option_t *option ="");
    virtual void Print(Option_t *option="") const;

protected:
    UInt_t m_status;        /// GAMMA filter status word
    UInt_t m_stage;         /// bit mask indicating which stages of the filter have been processed
    UInt_t m_energyValid;   /// Flag indicating the energy of the event was evaluated by the gamma filter
    Int_t  m_energyInLeus;  /// energy in units of LEUS (1/4 Mev)

private:

    ClassDef(LpaGammaFilter,1) // Gamma Filter
};
 
class LpaGammaFilterV0 : public LpaGammaFilter {
public:
    LpaGammaFilterV0() : LpaGammaFilter() { };
    LpaGammaFilterV0(UInt_t status) : LpaGammaFilter(status) { };
    virtual ~LpaGammaFilterV0() { }

private:

    ClassDef(LpaGammaFilterV0,1) // Gamma Filter Version 0
};

class LpaGammaFilterV1 : public LpaGammaFilter {
public:
    LpaGammaFilterV1() : LpaGammaFilter() { };
    LpaGammaFilterV1(UInt_t status) : LpaGammaFilter(status) { };
    virtual ~LpaGammaFilterV1() { }

private:

    ClassDef(LpaGammaFilterV1,1) // Gamma Filter Version 1 
};

class LpaGammaFilterV2 : public LpaGammaFilter {
public:
    LpaGammaFilterV2() : LpaGammaFilter() { };
    LpaGammaFilterV2(UInt_t status) : LpaGammaFilter(status) { };
    virtual ~LpaGammaFilterV2() { }

private:

    ClassDef(LpaGammaFilterV2,1) // Gamma Filter Version 2
};

class LpaHipFilter : public ILpaHandler
{
public:
    LpaHipFilter() : ILpaHandler(), m_status(0) {}
    LpaHipFilter(UInt_t status) : m_status(status) {}
    virtual ~LpaHipFilter() {}

    virtual void setStatusWord(UInt_t status) { m_status = status; }
    virtual UInt_t  getStatusWord() const { return m_status; };
    virtual UInt_t getStage() const { return 0; }
    
    virtual UInt_t getAllVetoBits() const;

    virtual void Clear(Option_t *option ="");
    virtual void Print(Option_t *option="") const;

protected:
    UInt_t  m_status;

private:

    ClassDef(LpaHipFilter,1) // Highly Ionizing Particle Filter output
};

class LpaHipFilterV0 : public LpaHipFilter {
public:
    LpaHipFilterV0() : LpaHipFilter() { };
    LpaHipFilterV0(UInt_t status) : LpaHipFilter(status) { };
    virtual ~LpaHipFilterV0() { }
    virtual UInt_t getStage() const {
       return m_status & LpaHipFilterV0::Stage_m;
    }

private:
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

    ClassDef(LpaHipFilterV0,1) // HIP Filter Version 0 
};

class LpaMipFilter : public ILpaHandler
{
public:
    LpaMipFilter(): ILpaHandler(), m_status(0) {}
    LpaMipFilter(UInt_t status) : m_status(status) {}
    virtual ~LpaMipFilter() {}

    virtual void setStatusWord(UInt_t status) { m_status = status; }
    virtual UInt_t  getStatusWord() const { return m_status; };
    virtual UInt_t getStage() const { return 0; }
    virtual UInt_t getAllVetoBits() const;


    virtual void Clear(Option_t *option ="");
    virtual void Print(Option_t *option="") const;

protected:
    UInt_t  m_status;

private:

    ClassDef(LpaMipFilter,1) // MIP filter output
};

class LpaMipFilterV0 : public LpaMipFilter {
public:
    LpaMipFilterV0() : LpaMipFilter() { };
    LpaMipFilterV0(UInt_t status) : LpaMipFilter(status) { };
    virtual ~LpaMipFilterV0() { }

    virtual UInt_t getStage() const { 
       return m_status & LpaMipFilterV0::Stage_m; 
    };

  private:
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


    ClassDef(LpaMipFilterV0,1) // MIP Filter Version 0 
};


class LpaDgnFilter : public ILpaHandler
{
public:
    LpaDgnFilter() : ILpaHandler(), m_status(0) {}
    LpaDgnFilter(UInt_t status) : m_status(status) {}
    virtual ~LpaDgnFilter() {}
 
    virtual void setStatusWord(UInt_t status) { m_status = status; }
    // If msb is set below then event is to be vetoed
    virtual UInt_t  getStatusWord() const { return m_status; };

    virtual UInt_t getStage() const { return 0; }

    virtual UInt_t getGemClasses() const { return 0; }

    virtual UInt_t getAllVetoBits() const;

    virtual void Clear(Option_t *option ="");
    virtual void Print(Option_t *option="") const;

protected:
    UInt_t  m_status;

private:
    ClassDef(LpaDgnFilter,1) // Diagnostic Filter output
};

class LpaDgnFilterV0 : public LpaDgnFilter {
public:
    LpaDgnFilterV0() : LpaDgnFilter() { };
    LpaDgnFilterV0(UInt_t status) : LpaDgnFilter(status) { };
    virtual ~LpaDgnFilterV0() { }
    virtual UInt_t getStage() const {
        return m_status & LpaDgnFilterV0::Stage_m;
    }

    virtual UInt_t getGemClasses() const {
        return (m_status & LpaDgnFilterV0::Gem_Classes_m) >> LpaDgnFilterV0::Gem_Classes_v;
    }

  private:
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

    ClassDef(LpaDgnFilterV0,1) // Diagnostic Filter output
};

  class LpaPassthruFilter : public ILpaHandler {
  public:
    LpaPassthruFilter() : ILpaHandler() { };

    LpaPassthruFilter(UInt_t status) : ILpaHandler() { m_status = status; }

    virtual ~LpaPassthruFilter() { }

    virtual void setStatusWord(UInt_t status) { m_status = status; }
 
    virtual UInt_t getStatusWord() const { return m_status; }

    virtual UInt_t getStage() const { return 0; }

    virtual void Clear(Option_t *option ="");
    virtual void Print(Option_t *option="") const;

  protected:

    UInt_t m_status;

  private:
    ClassDef(LpaPassthruFilter,1) // Passthru output

  };

 class LpaPassthruFilterV0 : public LpaPassthruFilter {
 public:
    LpaPassthruFilterV0() : LpaPassthruFilter() { };
    LpaPassthruFilterV0(UInt_t status) : LpaPassthruFilter(status) { };
    virtual ~LpaPassthruFilterV0() { }
    virtual UInt_t getStage() const { return 0; }

 private:
    ClassDef(LpaPassthruFilterV0,1) // Diagnostic Filter output
};

/** @class LpaFilterStatus
* @brief A container for the output of the OnboardFilter. Will contain a list of 
*        output objects for each of the filters run
*  
*/

class LpaHandler : public TObject  
{
public: 
    LpaHandler();
    
    virtual ~LpaHandler();
    
    /// clear lists, free pointers, etc., after read from / write to file
    void Clear(Option_t *option ="");
    void Fake( Int_t ievent, UInt_t rank, Float_t randNum ) ; // for tests
    Bool_t CompareInRange( const LpaHandler &, const std::string & name = "" ) const ; // for tests
    void Print(Option_t *option="") const;
    /// Assignment operator for convenience
    LpaHandler& operator=(const LpaHandler& rhs);
    
    void initialize();

    // Add results method
    void addHandler(const enums::Lsf::HandlerId &id,TObject* status);

    void addGamma(LpaGammaFilter* gam) { m_gamma = gam; }
    void addHip(LpaHipFilter* hip) { m_hip = hip;};
    void addMip(LpaMipFilter* mip) { m_mip = mip;};
    void addDgn(LpaDgnFilter* dgn) { m_dgn = dgn; };
    void addPassthru(LpaPassthruFilter* pass)
      { m_pass = pass; };
    void addAsc(ILpaHandler* asc) { m_asc = asc; }

    // Return results method
    const ILpaHandler* getHandler(const enums::Lsf::HandlerId &id) const;
    const LpaGammaFilter* getGammaFilter() const { return m_gamma; }
    const LpaHipFilter* getHipFilter() const { return m_hip; }
    const LpaMipFilter* getMipFilter() const { return m_mip; }
    const LpaDgnFilter* getDgnFilter() const { return m_dgn; }
    const LpaPassthruFilter* getPassthruFilter() const { return m_pass; }
    const ILpaHandler* getAsc() const { return m_asc; }

//    const TObjArray* getHandlerCol() const { return &m_lpaHandlerCol; }

private:
    /// An array of points defining the trajectory
  //  TObjArray m_lpaHandlerCol;
   LpaGammaFilter *m_gamma;
   LpaHipFilter *m_hip;
   LpaMipFilter *m_mip;
   LpaDgnFilter *m_dgn;
   LpaPassthruFilter *m_pass;
   ILpaHandler  *m_asc;
    
    ClassDef(LpaHandler,1) 
};


#endif

