//----------------------------------------------------------------------
//    
//    Implementation of the TkrFitHit Class
//
//      Original due to Jose Hernando-Angel circa 1997-1999
//      Re-written to combine both X and Y projections (2001) (Bill Atwood)
//      
//-----------------------------------------------------------------------

#include "Event/Recon/TkrRecon/TkrFitHit.h"

using namespace TkrRecon;

TkrFitHit TkrFitHit::changeType(TYPE typ)
{
    
    TkrFitHit hit;
    
    hit.m_type = typ;
    hit.m_par  = m_par;
    hit.m_cov  = m_cov;
    
    return hit;
}
