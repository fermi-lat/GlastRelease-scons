// $Header$
#ifndef CalIRFAlg_H
#define CalIRFAlg_H

// Include files
#include "GaudiKernel/Algorithm.h"
#include "CalGeometrySvc.h"
class IGlastDetSvc;


/*! \class CalIRFAlg
\brief Conversion calorimeter data from CsIData to CalRecLogs

  */

class CalIRFAlg : public Algorithm {

public:
  //! Constructor of this form must be provided
  CalIRFAlg(const std::string& name, ISvcLocator* pSvcLocator); 

  //! mandatory
  StatusCode initialize();
  //! mandatory
  StatusCode execute();
  //! mandatory
  StatusCode finalize();


private:
    // the GlastDetSvc used for access to detector info
    IGlastDetSvc*    m_detSvc;

	ICalGeometrySvc* m_CalGeo;

};


#endif // CalIRFAlg_H
