// $Header$
#include "CalibData/Cal/LightAtt.h"

namespace CalibData {

  void LightAtt::update(RangeBase* other)  {
    LightAtt* otherLightAtt = dynamic_cast<LightAtt* > (other);
    
    // check that otherLightAtt isn't 0 (dynamic cast worked)
    m_lightAtt = otherLightAtt->m_lightAtt;
    m_lightNorm = otherLightAtt->m_lightNorm;
  }
}
