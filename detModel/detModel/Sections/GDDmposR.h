#ifndef GDDMPOSR_H
#define GDDMPOSR_H
#include "detModel/Sections/GDDanyPosition.h"

/**
 * @author D.Favretto & R.Giannitrapani 
 */
class GDDmposR : GDDanyPosition {
public:

    GDDmposR(int pncopy,double pdR,positiontype ptypePosition);
    void setR0(double pr0);
    void setZPhi(double pZPhi[]);

private:    
    int ncopy;
    double r0; ///default 0
    double dR;
    double ZPhi[2]; ///default "0 0"
};
#endif //GDDMPOSR_H
