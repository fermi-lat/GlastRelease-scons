#ifndef calibRootData_ChargeScale_h
#define calibRootData_ChargeScale_h

#include "TObject.h"
#include "commonRootData/idents/TkrId.h"
#include "TClonesArray.h"

/** 
 *   @file ChargeScale.h
 *    defines Root classes used for charge scale calibration data
 *      ChargeScale Strip   is per-strip information
 *      ChargeScaleLayer   is per-Si layer information
 */

namespace calibRootData {

  /**
   *    @class ChargeScaleStrip
   *    @brief  Information associated with a single strip or gtfe
   */
  class ChargeScaleObj : public TObject {
  private:
    Int_t  m_Id;
    Float_t m_scale;
    Float_t m_error;
    Float_t m_chi2;
    Float_t m_df;

  public:
    ChargeScaleObj(Int_t id=-1, Float_t scale=0.0, Float_t error=0.0,
                     Float_t chi2 = 0.0, Float_t df = 0.0) :
      m_Id(id), m_scale(scale), m_error(error), m_chi2(chi2),
      m_df(df) {}

    UInt_t getId() const {return m_Id;}
    Float_t getScale() const {return m_scale;}
    Float_t getError() const {return m_error;}
    Float_t getChi2() const {return m_chi2;}
    Float_t getDf() const {return m_df;}

    void Clear(Option_t *option ="");
    void copy(const ChargeScaleObj& other);

    void init(UInt_t id, Float_t scale, Float_t error, Float_t chi2, 
              Float_t df);

  private:
    ClassDef(calibRootData::ChargeScaleObj, 1)
  };

  /**  @class ChargeScaleUnilayer
        Keeps track of one unilayer's worth of charge scale constants.
        Constants may be kept either per strip or per gtfe.  In creating
        this objects from an xml file, we don't know which it is until
        we see which kind of child object we've got, so don't allocate
        the TClonesArray until later, using resize(..)
  */
  class ChargeScaleUnilayer : public TObject {
  private:
    commonRootData::TkrId m_unilayerId; 
    UInt_t m_nChildren;
    TClonesArray* m_children;
    bool m_childIsStrip;
  public:
    ChargeScaleUnilayer(const commonRootData::TkrId& id);

    ChargeScaleUnilayer() : m_unilayerId(), m_nChildren(0), m_children(0),
                            m_childIsStrip(false) {};
    void setId(const commonRootData::TkrId& id) {m_unilayerId = id;}
    const commonRootData::TkrId& getId() const {return m_unilayerId;}
    TClonesArray* getChildren() {return m_children;}
    Int_t getNAlloc() const {return (m_children->LastIndex() + 1); }
    Int_t getNChildren() const {return m_nChildren; }
    bool  consistent() const {return (getNChildren() == getNAlloc());} 
    const ChargeScaleObj* getObj(Int_t n) const {
      return (ChargeScaleObj*)(m_children->At(n));
    }
    bool childIsStrip() const {return m_childIsStrip;}


    /** Make room for nChild worth of data.  Old data, if there is
     any, is lost.  We find out how many children at the same time
     we find out what kind they are
    */
    void resize(unsigned nChild, bool childIsStrip);

    void Clear(Option_t *option ="");

    /** Put strip or gtfe data in the proper place.
        @ret true if successful (strip or gtfe number within range); else false
    */
    Bool_t setChild(const ChargeScaleObj& child);
  private:
    ClassDef(calibRootData::ChargeScaleUnilayer, 1)
  };

}

#endif
