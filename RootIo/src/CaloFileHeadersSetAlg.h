
#include "GaudiKernel/Algorithm.h"

/** @class CaloFileHeadersSetAlg
 * @brief Prepare the calo specific attributes and set them in the relevant headers.
 *
 * @author David Chamont
 * $Header$
 */

class CaloFileHeadersSetAlg : public Algorithm  {
        
public:
    
    CaloFileHeadersSetAlg( const std::string& name, ISvcLocator* pSvcLocator ) ;
    
    /// 
    StatusCode initialize();
   
    /// 
    StatusCode execute();
    
    /// prepare the calo specific attributes and set them in the relevant headers
    StatusCode finalize();
    
} ;

