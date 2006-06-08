#ifndef calibRootData_Tot_h
#define calibRootData_Tot_h

#include "TObject.h"
#include "commonRootData/idents/TkrId.h"
#include "TClonesArray.h"

/** 
 *   @file Tot.h
 *    defines Root classes used for ToT calibration data
 *      ToTStrip   is per-strip information
 *      ToTLayer   is per-Si layer information
 */

namespace calibRootData {

  /**
   *    @class TotStrip
   *    @brief  Information associated with a single strip
   */
  class TotStrip : public TObject {
  private:
    Int_t  m_stripId;
    Float_t m_slope;
    Float_t m_intercept;
    Float_t m_quad;
    Float_t m_chi2;
    Float_t m_df;

  public:
    TotStrip(Int_t id=-1, Float_t slope=0.0, Float_t intercept=0.0,
      Float_t quad=0.0, Float_t chi2=0.0, Float_t df=0.0) : 
      m_stripId(id), m_slope(slope), m_intercept(intercept), m_quad(quad), 
      m_chi2(chi2), m_df(df) {}
    ~TotStrip() {}

    UInt_t getStripId() const {return m_stripId;}
    Float_t getSlope() const {return m_slope;}
    Float_t getIntercept() const {return m_intercept;}
    Float_t getQuad() const {return m_quad;}
    Float_t getChi2() const {return m_chi2;}
    Float_t getDf() const {return m_df;}

    void Clear(Option_t *option ="");
    void copy(const TotStrip& other);

    void init(UInt_t id, Float_t slope, Float_t intercept,
              Float_t quad, Float_t chi2, Float_t df);


  private:
    ClassDef(calibRootData::TotStrip, 1)
  };

  class TotUnilayer : public TObject {
  private:
    commonRootData::TkrId m_unilayerId; 
    UInt_t m_nStrips;
    TClonesArray* m_strips;
  public:
    TotUnilayer(const commonRootData::TkrId& id, UInt_t nStrips=1536);
    TotUnilayer() : m_unilayerId(), m_nStrips(0), m_strips(0) {};   
    virtual ~TotUnilayer();
    void setId(const commonRootData::TkrId& id) {m_unilayerId = id;}
    const commonRootData::TkrId& getId() const {return m_unilayerId;}
    TClonesArray* getStrips() {return m_strips;}
    Int_t getNAlloc() const {return (m_strips->LastIndex() + 1); }
    Int_t getNStrips() const {return m_nStrips; }
    bool  consistent() const {return (getNStrips() == getNAlloc());} 
    const TotStrip* getStrip(Int_t n) const {
      return (TotStrip*)(m_strips->At(n));
    }

    /** Make room for nStrips worth of data.  Old data, if there is
     any, is lost.
    */
    void resize(unsigned nStrips);

    void Clear(Option_t *option ="");

    /** Put strip data in the proper place.
        @ret true if successful (strip number within range); else false
    */
    Bool_t setStrip(const TotStrip& strip);
  private:
    ClassDef(calibRootData::TotUnilayer, 1)
  };

}

#endif
