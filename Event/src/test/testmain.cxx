#include "Event/TopLevel/EventModel.h"
#include "geometry/Vector.h"
#include "Event/Utilities/TimeStamp.h"
int main(){
    TimeStamp ghi(2009,12,19,11);
    std::cout << ghi << std::endl;
    return 0;
}