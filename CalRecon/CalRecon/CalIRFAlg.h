// $Header$
#ifndef CalIRFAlg_H
#define CalIRFAlg_H

// Include files
#include "Gaudi/Algorithm/Algorithm.h"

#include "reconstruction/SummaryData.h"
#include "CalRecon/GaudiGlastTuple.h"

// forward declarations
class IGlastDetSvc;
class CalRecon;
class GlastTuple;
namespace xml { class IFile; }


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

  StatusCode printNewNTuple();

private:
    // the GlastDetSvc used for access to detector info
    IGlastDetSvc*    m_detSvc;
    // ptr to the CalRecon object used to do the analysis
    CalRecon*    m_recon;

    // constants from the "instrument.xml" file.
    xml::IFile * m_ini;
    
    // sumamry object from glastsim creates a n-tuple
//    SummaryData<GlastTuple>* m_summary;
    SummaryData<GaudiGlastTuple>* m_gsummary;


};


#endif // CalIRFAlg_H
