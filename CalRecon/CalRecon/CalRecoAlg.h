// $Header$
#ifndef CalRecoAlg_H
#define CalRecoAlg_H

// Include files
#include "Gaudi/Algorithm/Algorithm.h"

#include "reconstruction/SummaryData.h"

// forward declarations
class IGlastDetSvc;
class CalRecon;
class GlastTuple;


/*! \class CalRecoAlg
\brief Calorimeter reconstuction

  */

class CalRecoAlg : public Algorithm {

public:
  //! Constructor of this form must be provided
  CalRecoAlg(const std::string& name, ISvcLocator* pSvcLocator); 

  //! mandatory
  StatusCode initialize();
  //! mandatory
  StatusCode execute();
  //! mandatory
  StatusCode finalize();

private:
    // the GlastDetSvc used for access to detector info
    IGlastDetSvc*    m_detSvc;
    // ptr to the CalRecon object used to do the analysis
    CalRecon*    m_recon;

    // sumamry object from glastsim creates a n-tuple
    SummaryData<GlastTuple>* m_summary; 
};


#endif // CalRecoAlg_H
