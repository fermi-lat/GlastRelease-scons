#include "detModel/Management/SectionsVisitor.h"
#include "detModel/Sections/Tube.h"
#include "detModel/Sections/BoundingBox.h"

namespace detModel{

  void Tube::buildBB(){
    getBBox()->setXDim(2*getRout());
    getBBox()->setYDim(2*getRout());
    getBBox()->setZDim(getZ());
  }
}
