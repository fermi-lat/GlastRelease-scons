#include "Event/TopLevel/EventModel.h"
#include "Event/Utilities/SkyDir.h"
#include "geometry/Vector.h"

int main(){
    Vector abc(0,0,1);
    SkyDir xyz(0,0,GALACTIC);
    //SkyDir xyz(abc);
    std::cout << " , ra =" << xyz.ra() << " , dec =" << xyz.dec() << std::endl;
    std::cout << " , l =" << xyz.l() << " , b =" << xyz.b() << std::endl;

    return 0;
}