// $Header$
#ifndef CalibData_TkrTot_h
#define CalibData_TkrTot_h

#include "CalibData/Tkr/TkrBase.h"
#include "CalibData/CalibModel.h"
#include "CalibData/Tkr/UniBase.h"
#include <vector>

/** 
       @file CalibData/TkrTot.h

       Defines structure of ToT data: array of 1536 strip objects
        per Si uniplane

       @author  J. Bogart
*/

class RootTkrTotCnv;

namespace CalibData {

  /**
     @class TkrTotStrip

     Just a pure data class.  
  */
  class TkrTotStrip {
  public:
  private:
    int   m_stripId;
    float m_slope;
    float m_intercept;
    float m_quad;
    float m_chi2;
    float m_df;

  public:
    TkrTotStrip(int id=-1, float slope=0.0, float intercept=0.0,
      float quad=0.0, float chi2=0.0, float df=0.0) : 
      m_stripId(id), m_slope(slope), m_intercept(intercept), m_quad(quad), 
      m_chi2(chi2), m_df(df) {}
    TkrTotStrip(const TkrTotStrip& other);
    void clear() {
      m_stripId = -1; 
      m_slope = m_intercept = m_quad = m_chi2 = m_df = 0.0;
    }
    int getStripId() const {return m_stripId;}
    float getSlope() const {return m_slope;}
    float getIntercept() const {return m_intercept;}
    float getQuad() const {return m_quad;}
    float getChi2() const {return m_chi2;}
    float getDf() const {return m_df;}

    void copy(const TkrTotStrip& other);

    void init(int id, float slope, float intercept,
              float quad, float chi2, float df);

  };

  /**
     @class TkrTotUni

     collection of strip objects belonging to the uniplane, plus
     enough information to identify the uniplane
  */
  class TkrTotUni : public UniBase {
  public:
    TkrTotUni(const idents::TkrId& id, int nStrips=0);
    TkrTotUni() : UniBase(), m_strips(0), m_nStrips(0) {}

    ~TkrTotUni() {if (m_strips) delete [] m_strips;};

    unsigned getNStrips() const {return m_nStrips;}
    const TkrTotStrip* getStripData(unsigned i) const {
      if (i > m_nStrips) return 0;
      return &m_strips[i];
    }

    // Is there a possibility that calibration files
    // will have "holes":  No data for stripId x but data for x+1 ? 
    void clearStrips();

    /// Resize array of strips, or just clear if sizes match
    void resize(unsigned n);
    /**
       Use strip id field to decide where to put strip. Success iff
       the id is in range
    */
    bool putStrip(const TkrTotStrip& strip);

    // reimplemented from UniBase
    virtual void update(UniBase* other);

  private:
    TkrTotStrip* m_strips;
    unsigned    m_nStrips;
  };



  /**
     class TkrTotCol

     TDS class for complete (charge-injection) calibration
  */
  class TkrTotCol : public TkrBase {
    friend class ::RootTkrTotCnv;

  public:
    TkrTotCol(unsigned nTowerRow=4, unsigned nTowerCol=4, unsigned nTray=19);

    /// Nothing left to do in destructor; it's all handled by base class
    ~TkrTotCol() {}

    /// The one routine of interest to applications: fetch constants for
    /// a particular strip
    const TkrTotStrip* getStripInfo(const idents::TkrId& id, unsigned iStrip) {
      TkrTotUni* uni = dynamic_cast<TkrTotUni*>(getUni(id));
      if (!uni) return 0;
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

  };     // end   TkrTotCol

}
#endif
