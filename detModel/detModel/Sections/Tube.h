#ifndef TUBE_H
#define TUBE_H
#include "detModel/Sections/Shape.h"

namespace detModel{

  class Tube :public Shape {
  private:
    double RIn;///default 0
    double ROut;///default 0
    double Z;///default 0
    double PhiBeg;///default 0
    double PhiEnd;///default 360
  };
}
#endif //TUBE_H
