
#include <RootConvert/Recon/CalClassParamsConvert.h>
#include <RootConvert/Utilities/Toolkit.h>

namespace RootPersistence {

void convert( const Event::CalClassParams & tdsClassParams, CalClassParams & rootClassParams ) {
    rootClassParams.init
      (  ) ;
 }
 
void convert( const CalClassParams & rootClassParams, Event::CalClassParams & tdsClassParams )
{
    tdsClassParams = Event::CalClassParams
      (  ) ;
}

}
