
#include "ICalEnergyCorr.h"
#include "ICalReconSvc.h"
#include "GaudiKernel/Algorithm.h"

// for implementation
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/DeclareFactoryEntries.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "Event/Recon/TkrRecon/TkrVertex.h"
#include "Event/Recon/CalRecon/CalCluster.h"
#include "Event/TopLevel/EventModel.h"

class IGlastDetSvc;

/**   
* @class CalEventEnergyAlg
*
* @brief Algorithm for reconstruction of energy and direction of incident particle
*
*
* Performs high level energy corrections
*
* The reconstruction here uses the CalClusterCol and produces
* instances of CalEnergyCorrResult.
* It tries to correct for energy leakage using two
* different methods:
*  -  Profile performs profile fitting
*  -  Leak performs a correction using the last layer correlation
* 
* See the description of Profile() and Leak() for details
* For most of the test beam energies Leak() should give better results.
* there can be a strong bias though because of miscalibration
*
* For a comparison one can see on the following plots the results of this
* method on R138 data and a MC run of 20 GeV positrons
* 
* \image html figurearticle2.gif
* \image latex figurearticle2.eps width=10cm
* \image html figuresim2.gif
* \image latex figuresim2.eps width=10cm
* \author Regis Terrier
* \author Alexandre Chekhtman
* \author Jose Hernando
* 
* \warning May not give sensible results when there is a large uncertainty in 
* the gains especially for the Leak method if there is a problem in the last 
* layer.
* \warning High energy corrections are intended for high energy 
* \todo Add low energy corrections 
*
* $Header$
*/


class CalEventEnergyAlg : public Algorithm
{
public:
    
    //! constructor
    CalEventEnergyAlg( const std::string & name, ISvcLocator * pSvcLocator ) ; 
    //! destructor
    virtual ~CalEventEnergyAlg() ; 
    
    virtual StatusCode initialize() ;

/*! Performs high energy corrections: see Profile() and Leak() for details
 */        
    StatusCode execute();

    StatusCode finalize() ;
        
private:
    
    //! correction tool names
    StringArrayProperty m_corrToolNames ;
    
    //! correction tools
    std::vector<ICalEnergyCorr*> m_corrTools ;
    
    //! Type of correction for primary corrected parameters
    unsigned int m_corrType;

    //! Set status bits depending on which iteration of algorithm
    unsigned int m_passBits;
    
    //! package service
    ICalReconSvc * m_calReconSvc ;
    
} ;


//==============================================
// IMPLEMENTATION
//==============================================


DECLARE_ALGORITHM_FACTORY(CalEventEnergyAlg) ;

using namespace Event;

CalEventEnergyAlg::CalEventEnergyAlg
 ( const std::string & name, ISvcLocator * pSvcLocator )
 : Algorithm(name,pSvcLocator), m_calReconSvc(0)
 {
    std::vector<std::string> corrToolNames ;
    unsigned int             corrType ;

    // Initial/default parameters depend on whether this would be the first pass
    // through CalRecon or the second. The first pass through should have the 
    // algorithm name "RawEnergy" assigned to it, the second (or subsequent) name
    // doesn't matter for initialization.
    // This branch if first pass
    if (name == "RawEnergy")
    {
        corrToolNames.push_back("CalRawEnergyCorr") ;
// FOR NEW TDS
//        corrType   = Event::CalCorToolResult::RAWENERGY ;
//        m_passBits = Event::CalEventEnergy::PASS_ONE ;
    }
    // This the default for the second iteration of CalRecon
    else
    {
        corrToolNames.push_back("CalLastLayerCorr");
        corrToolNames.push_back("CalProfileCorr");
        corrToolNames.push_back("CalValsCorr");
        corrToolNames.push_back("CalTkrLikelihoodCorr");
        corrToolNames.push_back("CalTransvOffsetCorr");
// FOR NEW TDS
//        corrType   = Event::CalCorToolResult::CALVALS ;
//        m_passBits = Event::CalEventEnergy::PASS_TWO ;
    }

    // Declare the properties with these defaults
    declareProperty("corrToolNames", m_corrToolNames = corrToolNames);
    declareProperty("finalCorrType", m_corrType      = corrType);
 }

CalEventEnergyAlg::~CalEventEnergyAlg()
 {
  m_corrTools.clear() ;
 }

StatusCode CalEventEnergyAlg::initialize()
 {
    MsgStream log(msgSvc(), name()) ;
    StatusCode sc = StatusCode::SUCCESS ;
    
    if (service("CalReconSvc",m_calReconSvc,true).isFailure())
     {
      log<<MSG::ERROR<<"Could not find CalReconSvc"<<endreq ;
      return StatusCode::FAILURE ;
     }

    const std::vector< std::string > & corrToolNames = m_corrToolNames ;
    std::vector< std::string >::const_iterator toolName ;
    ICalEnergyCorr * tool ;
    for ( toolName = corrToolNames.begin() ;
          toolName != corrToolNames.end() ;
          ++toolName ) {
        sc = toolSvc()->retrieveTool(*toolName,tool) ;
        if (sc.isFailure() ) {
            log<<MSG::ERROR<<" Unable to create "<<*toolName<<endreq ;
            return sc ;
        }
        else {
            m_corrTools.push_back(tool) ;
        }
    }

    return sc;
}

StatusCode CalEventEnergyAlg::execute()
{
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;
    
    log<<MSG::DEBUG<<"Begin"<<endreq ;
    m_calReconSvc->reviewEvent() ;
    
// FOR NEW TDS
//    // If no CalEnergy object (yet) then create one and register in TDS
//    Event::CalEventEnergy * calEnergy = SmartDataPtr<Event::CalEventEnergy>(eventSvc(),EventModel::CalRecon::CalEventEnergy);
//    if (calEnergy == 0)
//    {
//        calEnergy = new Event::CalEventEnergy();
//
//        if ((eventSvc()->registerObject(EventModel::CalRecon::CalEventEnergy, calEnergy)).isFailure())
//        {
//            log<<MSG::ERROR<<"Cannot register CalEventEnergy"<<endreq ;
//            return StatusCode::FAILURE ;
//        }
//    }

    // non fatal errors
    // if there's no clusters then CalEventEnergyAlg is doing nothing
  	if (!m_calReconSvc->getClusters())
      return StatusCode::SUCCESS ;
      
    // other errors are fatal
    if (m_calReconSvc->getStatus().isFailure())
      return StatusCode::FAILURE ;
      
    // loop over all correction tools
    int itool = 0 ;
    std::vector<ICalEnergyCorr *>::const_iterator tool ;
    for ( tool = m_corrTools.begin() ;
          tool != m_corrTools.end() ;
          ++tool, ++itool )
     {
        log<<MSG::DEBUG<<(*tool)->name()<<endreq ;
        if (((*tool)->doEnergyCorr()).isFailure())
         { 
          log<<MSG::ERROR<<"Failed"<<endreq ;
          sc = StatusCode::FAILURE ;
         }
     }
        
// FOR NEW TDS
//    // Set the pass number bits
//    calEnergy->setStatusBit(m_passBits);
//
//    // Go through and pick the "best" energy
//    // For now this consists of choosing the value from CalValsCorrTool
//    // It is envisaged that this will be replaced by a call to a tool/method which examines all 
//    // possibilites and chooses the best for the given event
//    for(Event::CalCorToolResultCol::iterator corIter = calEnergy->begin(); corIter != calEnergy->end(); corIter++)
//    {
//        Event::CalCorToolResult* corResult = *corIter;
//
//        if ((corResult->getStatusBits() & m_corrType) == m_corrType)
//        {
//            // Set the main energy parameters to the selected correction tool values
//            // Should also be setting a bit in CalEventEnergy to describe this?
//            calEnergy->setParams(corResult->getParams());
//            break;
//        }
//    }

    // print the reconstruction results for debugging
    m_calReconSvc->getClusters()->writeOut(log<<MSG::DEBUG);     
    
    log<<MSG::DEBUG<<"End"<<endreq ;
    return sc;
}

StatusCode CalEventEnergyAlg::finalize()
 { return StatusCode::SUCCESS ; }




