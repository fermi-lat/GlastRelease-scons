// $Header$
#ifndef CalibData_TkrTot_h
#define CalibData_TkrTot_h

#include "CalibData/Tkr/TkrBase.h"
#include "CalibData/CalibModel.h"
#include "CalibData/RangeBase.h"
#include <vector>

/** 
       @file CalibData/TkrTot.h

       Defines structure of ToT data: array of 1536 strip objects
        per Si uniplane

       @author  J. Bogart
*/

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
  class TkrTotUni : public RangeBase {
  public:
    TkrTotUni(const idents::TkrId& id, int nStrips=1536);
    TkrTotUni() : m_id(), m_strips(0), m_nStrips(0) {}

    ~TkrTotUni() {if (m_strips) delete [] m_strips;};

    //  Get functions
    unsigned getNStrips() const {return m_nStrips;}
    const TkrTotStrip* getStripData(unsigned i) const {
      if (i > m_nStrips) return 0;
      return &m_strips[i];
    }

    // reimplemented from RangeBase
    virtual void update(RangeBase* other);
    virtual void makeNew(RangeBase** ppNew);


  private:
    idents::TkrId m_id;
    TkrTotStrip* m_strips;
    unsigned    m_nStrips;
  };



  /**
     class TkrTotCol
  */
  class TkrTotCol : public TkrBase {
    friend class XmlTkrTotCnv;

  public:
    TkrTotCol(unsigned nTowerRow, unsigned nTowerCol, 
              unsigned nTray) : TkrBase(nTowerRow, nTowerCol, nTray) { }

    ~TkrTotCol();

    /// The one routine of interest to applications: fetch constants for
    /// a particular strip
    const TkrTotStrip* getStripInfo(const idents::TkrId& id, unsigned iStrip) {
      TkrTotUni* uni = dynamic_cast<TkrTotUni*>(getChannel(id));
      if (!uni) return 0;
      return uni->getStripData(iStrip);
    }

    // Reimplemented from DataObject
    virtual const CLID& clID() const {return classID(); }
    static const CLID&  classID();

    
    /// Reimplemented from TkrBase. Copy in TkrTotUni data
    virtual bool putChannel(RangeBase* data, const idents::TkrId& id, 
                            unsigned feChip=0);

    /// Reimplemented from TkrBase. Copy in TkrTotUni data
    virtual bool putChannel(RangeBase* data, unsigned towerRow, 
                            unsigned towerCol, unsigned tray, 
                            bool top, unsigned feChip=0) {
      return putChannel(data, idents::TkrId(towerCol, towerRow,
                                            tray, top) );
    }

    // Reimplemented from TkrBase. Copy in TkrTotUni data
    // or maybe we don't need to
    //    virtual StatusCode update(TkrBase& other, MsgStream* log);

    
  private:
    /// Pointers to (large) objects for uniplane-worth of data, allocated
    /// as needed during conversion from ROOT file.
    std::vector<TkrTotUni*> m_unis;

  };     // end   TkrTotCol

}
#endif
