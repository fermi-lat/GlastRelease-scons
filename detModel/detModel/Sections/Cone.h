#ifndef CONE_H
#define CONE_H
#include "detModel/Sections/Shape.h"

namespace detModel{

  class Cone :public Shape {
  private:
    double RIn1;///default 0
    double ROut1;///default 0
    double RIn2;///default 0
    double ROut2;///default 0
    double Z;///default 0
    double PhiBeg;///default 0
    double PhiEnd;///default 360
  };
}
#endif //CONE_H
