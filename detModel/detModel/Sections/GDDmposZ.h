#ifndef GDDMPOSZ_H
#define GDDMPOSZ_H
#include "detModel/Sections/GDDanyPosition.h"

/**
 * @author D.Favretto & R.Giannitrapani 
 */
class GDDmposZ : GDDanyPosition {
public:

    GDDmposZ(int pncopy,double pdz,positiontype ptypePosition);

    void setZ0(double pz0);
    void setRPhi(double pRPhi[]);
    void setXY(double pxy[]);

private:    
    int ncopy;
    double z0; ///default 0
    double dz;
    double RPhi[2]; ///default "0 0"
    double xy[2]; ///default "0 0"
}
;
#endif //GDDMPOSZ_H
