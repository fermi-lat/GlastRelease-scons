#include "Event/TopLevel/EventModel.h"
#include "Event/Utilities/SkyDir.h"
#include "geometry/Vector.h"
#include "Event/Utilities/TimeStamp.h"
int main(){
    SkyDir xyz(264.348,54.539,CELESTIAL);
    std::cout << " , ra =" << xyz.ra() << " , dec =" << xyz.dec() << std::endl;
    std::cout << " , l =" << xyz.l() << " , b =" << xyz.b() << std::endl;
    std::cout <<  "X =" << xyz.r().x() << " , Y =" << xyz.r().y() << " , Z =" << xyz.r().z() << std::endl;
    xyz().rotateX(1.);
    std::cout << " , ra =" << xyz.ra() << " , dec =" << xyz.dec() << std::endl;
    std::cout << " , l =" << xyz.l() << " , b =" << xyz.b() << std::endl;

    TimeStamp ghi(2009,12,19,11);
    std::cout << ghi << std::endl;
    return 0;
}