// $Header$
#ifndef FlightTest_h
#define FlightTest_h

// Include files
// Gaudi system includes
#include "GaudiKernel/Algorithm.h"

// ntuple interface
#include "ntupleWriterSvc/INTupleWriterSvc.h"


//------------------------------------------------------------------------------
/** 
A test algorithm, to allow devleopment, and to serve as a model for the "real" one.
  
*/
class FlightTest : public Algorithm {
public:
    FlightTest(const std::string& name, ISvcLocator* pSvcLocator);

    StatusCode initialize();
    StatusCode execute();
    StatusCode finalize();
    
private: 
    //! the GlastDetSvc used for access to detector info
    IGlastDetSvc*    m_detSvc; 

    //! constants from the "instrument.xml" file.
    xml::IFile * m_ini; 
    
    //! access the ntupleWriter service to write out to ROOT ntuples
    INTupleWriterSvc *m_ntupleWriteSvc;
    
    //! parameter to store the logical name of the ROOT file to write to
    std::string m_tupleName;
};
//------------------------------------------------------------------------
#endif
