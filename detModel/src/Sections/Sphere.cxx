#include "detModel/Management/SectionsVisitor.h"
#include "detModel/Sections/Sphere.h"
#include "detModel/Sections/BoundingBox.h"

namespace detModel{

  void Sphere::buildBB(){
    getBBox()->setXDim(2*m_rout);
    getBBox()->setYDim(2*m_rout);
    getBBox()->setZDim(2*m_rout);
  }
}
