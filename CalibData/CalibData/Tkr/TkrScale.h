// $Header$
#ifndef CalibData_TkrScale_h
#define CalibData_TkrScale_h

#include "CalibData/Tkr/TkrBase.h"
#include "CalibData/CalibModel.h"
#include "CalibData/Tkr/UniBase.h"
#include <vector>

/** 
       @file CalibData/TkrScale.h

       Defines structure of charge scale data: array of 1536 (if by strip)
       or 24 (if by gtfe) objects.  However interface in either case
       is to return per-strip constants.

       @author  J. Bogart
*/

class RootTkrScaleCnv;

namespace CalibData
 {

  /**
      @class TkrScaleObj 
      A pure data class
  **/

  class TkrScaleObj {
  public:
  private:
    int   m_id;
    float m_scale;
    float m_error;
    float m_chi2;
    float m_df;

  public:
    TkrScaleObj(int id=-1, float scale=0.0, float error=0.0,
      float chi2=0.0,  float df=0.0) : 
      m_id(id), m_scale(scale), m_error(error), m_chi2(chi2), m_df(df) {}
    TkrScaleObj(const TkrScaleObj& other);
    void clear() {
      m_id = -1; 
      m_scale = m_error = m_chi2 = m_df = 0.0;
    }
    int getId() const {return m_id;}
    float getScale() const {return m_scale;}
    float getError() const {return m_error;}
    float getChi2() const {return m_chi2;}
    float getDf() const {return m_df;}

    void copy(const TkrScaleObj& other);

    void init(int id, float scale, float error, float chi2, float df);

  };


  /**
     @class TkrScaleUni

     collection of strip objects belonging to the uniplane, plus
     enough information to identify the uniplane
  */
  class TkrScaleUni : public UniBase {
  public:
    TkrScaleUni(const idents::TkrId& id, int nObjs=0, bool perStrip = false);
    TkrScaleUni() : UniBase(), m_objs(0), m_nObjs(0), m_perStrip(false) {}

    ~TkrScaleUni() {if (m_objs) delete [] m_objs;};

    unsigned getNStrips() const {return m_nStrips;}
    unsigned getNObjs() const {return m_nObjs;}
    bool     perStrip() const {return m_perStrip;}

    /** Callers note: object returned is passed by value.  Since "real"
        information may be kept only per gtfe, object for strip must
         be manufactured
    */
    TkrScaleObj getStripData(unsigned i) const;

    // Is there a possibility that calibration files
    // will have "holes":  No data for stripId x but data for x+1 ? 
    void clearObjs();  

    /// Resize array of objects, or just clear if sizes match
    void resize(unsigned n, bool perStrip);
    /**
       Use strip id field to decide where to put strip. Success iff
       the id is in range
    */
    bool putObj(const TkrScaleObj& obj);

    // reimplemented from UniBase
    virtual void update(UniBase* other);

  private:
    TkrScaleObj*  m_objs;
    unsigned    m_nStrips;
    unsigned    m_nObjs;
    bool        m_perStrip;
  };


  /**
     class TkrScaleCol

     TDS class for complete (charge-injection) calibration
  */
  class TkrScaleCol : public TkrBase {
    friend class ::RootTkrScaleCnv;

  public:
    TkrScaleCol(unsigned nTowerRow=4, unsigned nTowerCol=4, unsigned nTray=19);

    /// Nothing left to do in destructor; it's all handled by base class
    ~TkrScaleCol() {}

    /// The one routine of interest to applications: fetch constants for
    /// a particular strip
    TkrScaleObj getStripInfo(const idents::TkrId& id, unsigned iStrip) {
      TkrScaleUni* uni = dynamic_cast<TkrScaleUni*>(getUni(id));
      if (!uni) return TkrScaleObj();
      return uni->getStripData(iStrip);
    }

    // Reimplemented from DataObject
    virtual const CLID& clID() const {return classID(); }
    static const CLID&  classID();

    // Reimplemented from TkrBase. All we need to do is check that
    // pointer is really of right type, then invoke base implementation.
    virtual bool putUni(UniBase* data, const idents::TkrId& id);

  private:
    // no new private data over what's in TkrBase

  };     // end   TkrScaleCol

}

#endif
