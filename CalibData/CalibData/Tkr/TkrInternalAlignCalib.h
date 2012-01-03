// $Header$
#ifndef CalibData_TkrInternalAlignCalib_h
#define CalibData_TkrInternalAlignCalib_h

#include "CalibData/CalibModel.h"
#include "CalibData/RangeBase.h"
#include "CalibData/CalibBase.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "GaudiKernel/StatusCode.h"
#include <vector>
#include <map>

class XmlTkrInternalAlignCnv;

namespace CalibData {

  // Or should this just be nested in TkrInternalAlignCalib ?
  class TkrInternalAlign : public RangeBase {
  public: 
    TkrInternalAlign();

    TkrInternalAlign(const CLHEP::Hep3Vector& disp, 
                     const CLHEP::Hep3Vector& rot)
    { putAlign(disp, rot); }
    virtual ~TkrInternalAlign() {}

    void getAlign(CLHEP::Hep3Vector& disp, CLHEP::Hep3Vector& rot);
    void putAlign(const CLHEP::Hep3Vector& disp, const CLHEP::Hep3Vector& rot);

  private:
    CLHEP::Hep3Vector m_disp;
    CLHEP::Hep3Vector m_rot;
                  
  };

  class TkrInternalAlignCalib : public CalibBase {
    friend class ::XmlTkrInternalAlignCnv;         // to be written

  public:
    TkrInternalAlignCalib()    {    };
    virtual ~TkrInternalAlignCalib();

    StatusCode getTrayAlign(unsigned towerId, unsigned trayId,
                            CLHEP::Hep3Vector& disp, CLHEP::Hep3Vector& rot);

    StatusCode getFaceAlign(unsigned towerId, unsigned trayId, unsigned faceId,
                            CLHEP::Hep3Vector& disp, CLHEP::Hep3Vector& rot);

    StatusCode getLadderAlign(unsigned towerId, unsigned trayId, 
                              unsigned faceId, unsigned ladderId,
                              CLHEP::Hep3Vector& disp, CLHEP::Hep3Vector& rot);

    StatusCode getWaferAlign(unsigned towerId, unsigned trayId,unsigned faceId,
                             unsigned ladderId, unsigned waferId,
                             CLHEP::Hep3Vector& disp, CLHEP::Hep3Vector& rot);

    // Re-implemented from CalibBase
    virtual StatusCode update(CalibBase& other, MsgStream* log);
    static const CLID&  classID();
                      
  private:
    unsigned makeTowerKey(unsigned tower) const;
    unsigned makeTrayKey(unsigned tower, unsigned tray) const;
    unsigned makeFaceKey(unsigned tower, unsigned tray, unsigned face) const;
    unsigned makeLadderKey(unsigned tower, unsigned tray, unsigned face,
                           unsigned ladder) const;
    unsigned makeWaferKey(unsigned tower, unsigned tray, unsigned face,
                          unsigned ladder, unsigned wafer) const;

    /**
       Given key (which encodes item information, get associated constants)
     */
    StatusCode getAlign(unsigned key, CLHEP::Hep3Vector& disp,
                        CLHEP::Hep3Vector& rot);

    /**
       Given key (which encodes item information).  Make new item if need
       be; else replace existing.
     */
    StatusCode putAlign(unsigned key, const CLHEP::Hep3Vector& disp,
                        const CLHEP::Hep3Vector& rot);

    void clean();

    /**
     Check that proposed item will have a key compatible with bit
     allocations.  
     @return  true means ok (within limits)
    */ 
    bool checkLimits(unsigned* tower, unsigned* tray=0, unsigned* face=0,
                     unsigned* ladder=0, unsigned* wafer=0);
    // Form key similar to TkrId; this is the key used to find disp, rot

    typedef std::map<unsigned int, TkrInternalAlign*> AlignMap;
    AlignMap m_items;
    //    std::map<unsigned int, TkrInternalAlign*> m_items
  };
}
#endif
