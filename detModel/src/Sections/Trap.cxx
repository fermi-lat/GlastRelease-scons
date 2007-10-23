#include "detModel/Management/SectionsVisitor.h"
#include "detModel/Sections/Trap.h"
#include "detModel/Sections/BoundingBox.h"

namespace detModel{

  void Trap::Accept(SectionsVisitor* v){
    v->visitTrap(this);
  }


  void Trap::buildBB(){

    double minX, maxX;
    minX = (-m_x1/2 < (-m_x2/2 + m_xDiff)) ? -m_x1/2 : -m_x2/2 + m_xDiff;
    maxX = (m_x1/2 > (m_x2/2 + m_xDiff) ) ? m_x1/2 : m_x2/2 + m_xDiff;
            getBBox()->setXDim(maxX - minX);
    getBBox()->setYDim(m_y);
    getBBox()->setZDim(m_z);
  }
  
}
