// $Header$
#ifndef CalRecoAlg_H
#define CalRecoAlg_H

// Include files
#include "Gaudi/Algorithm/Algorithm.h"

// forward declarations
class IGlastDetSvc;

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
};


#endif // CalRecoAlg_H
