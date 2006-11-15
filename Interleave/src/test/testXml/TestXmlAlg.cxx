/** @file TestXmlAlg.cxx
    @brief Used for test program

 $Header$

*/

// Include files
#include "../../XmlFetchEvents.h"

//#include "../BackgroundSelection.h"
// Gaudi system includes
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/SmartRefVector.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/Algorithm.h"

#include "TChain.h"

#include <string>

/*! \class TestXmlAlg
\brief 

*/

class TestXmlAlg : public Algorithm {

public:
    //! Constructor of this form must be provided
    TestXmlAlg(const std::string& name, ISvcLocator* pSvcLocator); 

    StatusCode initialize();
    StatusCode execute();
    StatusCode finalize();


private:
 
    XmlFetchEvents *m_fetch;
    std::string m_xmlFile;
};


static const AlgFactory<TestXmlAlg>  Factory;
//const IAlgFactory& CRTestAlgFactory = Factory;

//------------------------------------------------------------------------------
//
TestXmlAlg::TestXmlAlg(const std::string& name, ISvcLocator* pSvcLocator) 
: Algorithm(name, pSvcLocator)
{
    declareProperty("xmlFile",m_xmlFile="");

}

//------------------------------------------------------------------------------
/*! */
StatusCode TestXmlAlg::initialize() {

    StatusCode sc = StatusCode::SUCCESS;

    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "initialize" << endreq;

    // Use the Job options service to set the Algorithm's parameters
    setProperties();

    m_fetch = new XmlFetchEvents(m_xmlFile, "PtMagLat");

    return sc;
}


//------------------------------------------------------------------------------
StatusCode TestXmlAlg::execute() {

    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream   log( msgSvc(), name() );    

    double triggerRate = m_fetch->getAttributeValue("triggerRate", 5.0);
    double downlinkRate = m_fetch->getAttributeValue("downlinkRate",5.0);
    if ( (triggerRate == m_fetch->m_badVal) || (downlinkRate == m_fetch->m_badVal) )
        log << MSG::WARNING << "Failed to retrieve triggerRate or downlinkRate" << endreq;
    TChain *ch = new TChain();  // not specifying a TTree name, as our input tree names will come from XML
    int status = m_fetch->getFiles(5.0, ch);
    if (status == 1)
        log << MSG::WARNING << "At least one file failed to be added to the TChain" << endreq;
    else if (status == -1)
        log << MSG::WARNING << "No files were found" << endreq;
    
    return sc;
}


//------------------------------------------------------------------------------
StatusCode TestXmlAlg::finalize() {

    if (m_fetch) delete m_fetch;
    m_fetch = 0;
    return StatusCode::SUCCESS;
}






