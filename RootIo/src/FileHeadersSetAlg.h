
#include "GaudiKernel/Algorithm.h"

/** @class FileHeadersSetAlg
 * @brief Prepare the common attributes and set them in available headers.
 *
 * @author David Chamont
 * $Header$
 */

class FileHeadersSetAlg : public Algorithm {	

public:

    FileHeadersSetAlg( const std::string& name, ISvcLocator* pSvcLocator ) ;
    
    ///
    StatusCode initialize();
   
    ///
    StatusCode execute();
    
    /// prepare the common attributes and set them in available headers
    StatusCode finalize();
    
} ;

