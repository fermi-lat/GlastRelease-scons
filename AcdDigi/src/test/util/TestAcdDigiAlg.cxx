#define AcdDigi_TestAcdDigiAlg_CXX

// File and Version Information
// $Header$
// Description:
// Test for AcdUtil class. 

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"

#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "Event/TopLevel/EventModel.h"
#include "Event/TopLevel/Event.h"
#include "GaudiKernel/ObjectVector.h"

#include "GaudiKernel/Algorithm.h"

#include "../../AcdDigiUtil.h"

/** @class TestAcdDigiAlg.cpp
 * @brief AcdDigiUtil test algorithm
 *
 * Exercise all of AcdDigiUtil to be sure that the methods function properly.
 *
 * $Header$
 */

class TestAcdDigiAlg : public Algorithm {

public:
  TestAcdDigiAlg(const std::string& name, ISvcLocator* pSvcLocator); 

  StatusCode initialize();

  StatusCode execute();

  StatusCode finalize();

private:

    /// local member for access to AcdDigiUtil
    AcdDigiUtil m_util;

    /// input XML file containing parameters for Digitization
    std::string	m_xmlFile;

    /// vectors to store results of sampling
    std::vector<long> m_poisson;
    std::vector<float> m_gauss;
};

static const AlgFactory<TestAcdDigiAlg>  Factory;
const IAlgFactory& TestAcdDigiAlgFactory = Factory;

TestAcdDigiAlg::TestAcdDigiAlg(const std::string& name, ISvcLocator* pSvcLocator) :
Algorithm(name, pSvcLocator) {
    
    // Declare the properties that may be set in the job options file
    declareProperty ("xmlFile", m_xmlFile="$(ACDDIGIROOT)/xml/acdDigi.xml");
}

StatusCode TestAcdDigiAlg::initialize() {
    
    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "initialize" << endreq;
    StatusCode sc = StatusCode::SUCCESS;
    
    // Use the Job options service to set the Algorithm's parameters
    setProperties();

    m_util.getParameters(m_xmlFile);

    return sc;
}


StatusCode TestAcdDigiAlg::execute() {
    // Purpose and Method:  Exercise the methods contained in AcdDigiUtil
    //  When sampling the poisson and gaussian distributions, store the results to
    //  the appropriate vector.

    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream   log( msgSvc(), name() );

    float mips = m_util.convertMevToMips(1.0);
    log << MSG::DEBUG << "1.0 MeV = " << mips << " MIPs" << endreq;
    float pmtA_mips = mips;
    mips = m_util.convertMevToMips(20.0);
    log <<  MSG::DEBUG << "20.0 MeV = " << mips << " MIPs" << endreq;
    float pmtB_mips = mips;

    idents::AcdId id(0, 1, 1, 1);

    unsigned int pmtA_pe, pmtB_pe;

    m_util.convertMipsToPhotoElectrons(id, pmtA_mips, pmtA_pe,
        pmtB_mips, pmtB_pe);

    log << MSG::DEBUG << pmtA_mips << " is " << pmtA_pe << " pes" << endreq;
    log << MSG::DEBUG << pmtB_mips << " is " << pmtB_pe << " pes" << endreq;

    float pmtA2_mips, pmtB2_mips;

    m_util.convertPhotoElectronsToMips(id, pmtA_pe, pmtA2_mips,
        pmtB_pe, pmtB2_mips);

    log << MSG::DEBUG << pmtA_pe << " is " << pmtA2_mips << " MIPs" << endreq;
    log << MSG::DEBUG << pmtB_pe << " is " << pmtB2_mips << " MIPs" << endreq;

    float mipsToFullScaleA, mipsToFullScaleB;

    m_util.calcMipsToFullScale(id, pmtA_mips, pmtA_pe, 
        mipsToFullScaleA, pmtB_mips, pmtB_pe, mipsToFullScaleB );

    log << MSG::DEBUG << "mipsToFullScaleA " << mipsToFullScaleA << endreq;
    log << MSG::DEBUG << "mipsToFullScaleB " << mipsToFullScaleB << endreq;

    unsigned short pmtA_pha = m_util.convertMipsToPha(pmtA_mips, mipsToFullScaleA);
    unsigned short pmtB_pha = m_util.convertMipsToPha(pmtB_mips, mipsToFullScaleB);

    log << MSG::DEBUG << "pmtA_pha " << pmtA_pha << endreq;
    log << MSG::DEBUG << "pmtB_pha " << pmtB_pha << endreq;

    m_poisson.push_back(m_util.shootPoisson(5));

    m_gauss.push_back(m_util.shootGaussian(1.0));

    return sc;
}


StatusCode TestAcdDigiAlg::finalize() {
    // Purpose and Method:  Check the results of sampling the Poisson and Gaussian
    //   distributions.  Calculate the mean and standard deviation.
    
    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "finalize" << endreq;
    
    float poissonMean = 0.;
    std::vector<long>::iterator poissonIt;
    for (poissonIt = m_poisson.begin(); poissonIt != m_poisson.end(); poissonIt++) {
        poissonMean += (*poissonIt);
    }
    float poissonStdDev = 0.0;
    poissonMean /= m_poisson.size();
    for (poissonIt = m_poisson.begin(); poissonIt != m_poisson.end(); poissonIt++) {
        poissonStdDev += (*poissonIt - poissonMean) * (*poissonIt - poissonMean);
    }
    poissonStdDev /= m_poisson.size();
    log << MSG::DEBUG << "Poisson with actual mean 5 returned mean: " << poissonMean
        << " stdDev: " << poissonStdDev << endreq;

    std::vector<float>::iterator gaussIt;
    float gaussMean = 0.;
    for (gaussIt = m_gauss.begin(); gaussIt != m_gauss.end(); gaussIt++) {
        gaussMean += (*gaussIt);
    }
    float gaussStdDev = 0.0;
    gaussMean /= m_gauss.size();
    for (gaussIt = m_gauss.begin(); gaussIt != m_gauss.end(); gaussIt++) {
        gaussStdDev += (*gaussIt - gaussMean) * (*gaussIt - gaussMean);
    }
    gaussStdDev /= m_gauss.size();

    log << MSG::DEBUG << "Gaussian with actual stdDev 1.0 returned mean: " << gaussMean
        << " stdDev: " << gaussStdDev << endreq;


    return StatusCode::SUCCESS;
}






