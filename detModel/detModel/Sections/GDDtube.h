#ifndef GDDTUBE_H
#define GDDTUBE_H
#include "detModel/Sections/GDDshape.h"
class GDDtube :public GDDshape {
    private:
        double RIn;///default 0
        double ROut;///default 0
        double Z;///default 0
        double PhiBeg;///default 0
        double PhiEnd;///default 360
};
#endif //GDDTUBE_H
