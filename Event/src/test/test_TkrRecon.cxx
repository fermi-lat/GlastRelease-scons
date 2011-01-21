#include "Event/Recon/TkrRecon/TkrTrackParams.h"
#include "Event/Recon/TkrRecon/TkrTrackHit.h"
#include "GaudiKernel/ContainedObject.h"
#include "GaudiKernel/ObjectList.h"
#include <vector>
#include <iostream>
#include <string>



int main() 
{
    // Not really sure what a test program would do now, these classes only hold data...
    Event::TkrTrackParams params(1.,0.5,1.,0.5,0.100,0.,0.,0.,0.01,0.,0.,0.100,0.,0.01);
    Event::TkrTrackHit    hit;

    hit.getTrackParams(Event::TkrTrackHit::MEASURED) = params;

    std::cout << "Testing TkrTrackHit:  " << std::endl;
    std::cout << "  Measured x position:" << hit.getTrackParams(Event::TkrTrackHit::MEASURED)(1) << std::endl;
    std::cout << "  Measured x slope:   " << hit.getTrackParams(Event::TkrTrackHit::MEASURED)(2) << std::endl;
    std::cout << "  Measured y position:" << hit.getTrackParams(Event::TkrTrackHit::MEASURED)(3) << std::endl;
    std::cout << "  Measured y slope:   " << hit.getTrackParams(Event::TkrTrackHit::MEASURED)(4) << std::endl;
    std::cout << "  Measured x error:   " << hit.getTrackParams(Event::TkrTrackHit::MEASURED)(1,1) << std::endl;

        return 0;
}
