#include "Event/TopLevel/EventModel.h"
#include "Event/Utilities/SkyDir.h"
#include "geometry/Vector.h"
#include "Event/Utilities/TimeStamp.h"
int main(){
    //Vector abc(0,0,1);
    SkyDir xyz(0,0,CELESTIAL);
    //SkyDir xyz2(abc);
    std::cout << " , ra =" << xyz.ra() << " , dec =" << xyz.dec() << std::endl;
    std::cout << " , l =" << xyz.l() << " , b =" << xyz.b() << std::endl;
    xyz().rotateX(1.);
    std::cout << " , ra =" << xyz.ra() << " , dec =" << xyz.dec() << std::endl;
    std::cout << " , l =" << xyz.l() << " , b =" << xyz.b() << std::endl;

    TimeStamp ghi(2009,12,19,11);
    std::cout << ghi << std::endl;
    return 0;
}