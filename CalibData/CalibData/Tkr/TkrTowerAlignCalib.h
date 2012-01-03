// $Header$
#ifndef CalibData_TkrTowerAlignCalib_h
#define CalibData_TkrTowerAlignCalib_h

#include "CalibData/CalibModel.h"
#include "CalibData/RangeBase.h"
#include "CalibData/CalibBase.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "GaudiKernel/StatusCode.h"
#include <vector>

class XmlTkrTowerAlignCnv;

namespace CalibData {

  class TkrTowerAlign : public RangeBase {
  public: 
    TkrTowerAlign(const CLHEP::Hep3Vector& disp=CLHEP::Hep3Vector(),
                  const CLHEP::Hep3Vector& rot=CLHEP::Hep3Vector()) :
      m_disp(disp), m_rot(rot) {};
    virtual ~TkrTowerAlign() {}
    void getAlign(CLHEP::Hep3Vector& disp, CLHEP::Hep3Vector& rot);
    void putAlign(const CLHEP::Hep3Vector& disp, const CLHEP::Hep3Vector& rot);

    virtual void update(RangeBase* other);
  private:
    CLHEP::Hep3Vector m_disp;
    CLHEP::Hep3Vector m_rot;
  };

  class TkrTowerAlignCalib : public CalibBase {
    friend class ::XmlTkrTowerAlignCnv;         // to be written

  public:
    TkrTowerAlignCalib(unsigned maxTowerId=0) {
      m_towers.resize(maxTowerId + 1);
    }
    virtual ~TkrTowerAlignCalib() {
      m_towers.clear();
    }

    StatusCode getTowerAlign(unsigned id, CLHEP::Hep3Vector& disp,
                             CLHEP::Hep3Vector& rot) {
      if (id > m_towers.size()) return StatusCode::FAILURE;
      m_towers[id].getAlign(disp, rot);
      return StatusCode::SUCCESS;
    }

    // Re-implemented from CalibBase
    virtual StatusCode update(CalibBase& other, MsgStream* log);
    static const CLID&  classID();


                      
  private:
    StatusCode putAlign(unsigned id, const CLHEP::Hep3Vector& disp,
                        const CLHEP::Hep3Vector& rot) {
      if (id > m_towers.size()) return StatusCode::FAILURE;
      m_towers[id].putAlign(disp, rot);
      return StatusCode::SUCCESS;
    }
    std::vector<TkrTowerAlign> m_towers;
  };
}
#endif
