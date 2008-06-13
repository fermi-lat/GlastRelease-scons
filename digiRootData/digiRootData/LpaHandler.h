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

    void set(UInt_t status, UInt_t stage, UInt_t energyValid, Int_t energyInLeus) {
        m_status = status;
        m_stage = stage;
        m_energyValid = energyValid;
        m_energyInLeus = energyInLeus;
    }

    virtual UInt_t getStatusWord() const { return m_status; };
    UInt_t getStage() const { return m_stage; };
    UInt_t getEnergyValid() const { return m_energyValid; }
    Int_t getEnergyInLeus() const { return m_energyInLeus; }

    UInt_t getAllVetoBits() const;

    void Clear(Option_t *option ="");
    void Print(Option_t *option="") const;

private:
    UInt_t m_status;        /// GAMMA filter status word
    UInt_t m_stage;         /// bit mask indicating which stages of the filter have been processed
    UInt_t m_energyValid;   /// Flag indicating the energy of the event was evaluated by the gamma filter
    Int_t  m_energyInLeus;  /// energy in units of LEUS (1/4 Mev)

    ClassDef(LpaGammaFilter,1) // Gamma Filter
};

class LpaHipFilter : public ILpaHandler
{
public:
    LpaHipFilter() : ILpaHandler(), m_status(0) {}
    LpaHipFilter(UInt_t status) : m_status(status) {}
    virtual ~LpaHipFilter() {}

    void setStatusWord(UInt_t status) { m_status = status; }
    virtual UInt_t  getStatusWord() const { return m_status; };
    
    UInt_t getAllVetoBits() const;

    void Clear(Option_t *option ="");
    void Print(Option_t *option="") const;

private:
    UInt_t  m_status;

    ClassDef(LpaHipFilter,1) // Highly Ionizing Particle Filter output
};

class LpaMipFilter : public ILpaHandler
{
public:
    LpaMipFilter(): ILpaHandler(), m_status(0) {}
    LpaMipFilter(UInt_t status) : m_status(status) {}
    virtual ~LpaMipFilter() {}

    void setStatusWord(UInt_t status) { m_status = status; }
    virtual UInt_t  getStatusWord() const { return m_status; };

    UInt_t getAllVetoBits() const;


    void Clear(Option_t *option ="");
    void Print(Option_t *option="") const;

private:
    UInt_t  m_status;

    ClassDef(LpaMipFilter,1) // MIP filter output
};

class LpaDgnFilter : public ILpaHandler
{
public:
    LpaDgnFilter() : ILpaHandler(), m_status(0) {}
    LpaDgnFilter(UInt_t status) : m_status(status) {}
    virtual ~LpaDgnFilter() {}
 
    void setStatusWord(UInt_t status) { m_status = status; }
    // If msb is set below then event is to be vetoed
    virtual UInt_t  getStatusWord() const { return m_status; };

    UInt_t getAllVetoBits() const;

    void Clear(Option_t *option ="");
    void Print(Option_t *option="") const;

private:
    UInt_t  m_status;

    ClassDef(LpaDgnFilter,1) // Diagnostic Filter output
};


  class LpaPassthruFilter : public ILpaHandler {
  public:
    LpaPassthruFilter() : ILpaHandler() { };

    LpaPassthruFilter(UInt_t status) : ILpaHandler() { m_status = status; }

    void setStatusWord(UInt_t status) { m_status = status; }
 
    virtual UInt_t getStatusWord() const { return m_status; }

    void Clear(Option_t *option ="");
    void Print(Option_t *option="") const;

  private:

    UInt_t m_status;

    ClassDef(LpaPassthruFilter,1) // Passthru output

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

