// Gaudi system includes
#include "GaudiKernel/Algorithm.h"

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/IAlgManager.h"

#include "GuiSvc/IGuiSvc.h"

#include "gui/SubMenu.h"
#include "gui/GuiMgr.h"
#include "gui/Command.h"
#include "gui/Menu.h"

// TDS class declarations: input data, and McParticle tree
#include <cassert>

/** @class Replay 
 * @brief Call (replay?) a given algorithm by name, settable from job options
 * 
 * $Header$
*/
class Replay : public Algorithm {
public:
    Replay(const std::string& name, ISvcLocator* pSvcLocator);

    StatusCode initialize();
    StatusCode execute();
    StatusCode finalize(){return StatusCode::SUCCESS;}
    
private:
    // Data Members
    gui::GuiMgr* m_guiMgr;

    IAlgorithm* m_alg;
    StringProperty m_algname;

};

static const AlgFactory<Replay>  Factory;
const IAlgFactory& ReplayFactory = Factory;

Replay::Replay(const std::string& name, ISvcLocator* pSvcLocator)
:Algorithm(name, pSvcLocator)
{
    // declare properties with setProperties calls
    declareProperty("algName",  m_algname="TkrReconAlg");
    
}

//! set parameters and attach to various perhaps useful services.
StatusCode Replay::initialize(){
    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "initialize" << endreq;
    
    // Use the Job options service to set the Algorithm's parameters
    setProperties();
    
       
   // find the given algorithm 
    IAlgManager* theAlgMgr;
    if ( (serviceLocator( )->getService( "ApplicationMgr",
        IID_IAlgManager,
        (IInterface*&)theAlgMgr )).isFailure() ) {
        log << MSG::ERROR << " no application manager???" << endreq;
        return StatusCode::FAILURE;
        
    }
    

    if( (theAlgMgr->getAlgorithm(m_algname,  m_alg)).isFailure() ) {
        log << MSG::WARNING 
            << "Did not find the Algorithm "<< m_algname  << endreq;
        return StatusCode::FAILURE;
    }

        // get the  Gui service
    //
    IGuiSvc* guiSvc=0;
    
    if ( service("GuiSvc", guiSvc).isFailure() ){
        log << MSG::ERROR << "No GuiSvc, no display" << endreq;
        return StatusCode::SUCCESS; // go on anyway     
    }
    m_guiMgr = guiSvc->guiMgr();
    
    
    gui::GuiMgr* guiMgr=guiSvc->guiMgr() ;
    gui::Menu& menu = guiMgr->menu();

    // define the recall class right here: all it does is call the alg.
    class ReplayCommand: public gui::Command {
    public:
        ReplayCommand(IAlgorithm* alg):m_alg(alg){}
        void execute(){m_alg->execute();}
    private:
        IAlgorithm * m_alg;
    };
    menu.add_button(m_algname, new ReplayCommand(m_alg));

   return sc;
}

StatusCode Replay::execute()
{
    
    return StatusCode::SUCCESS; //m_alg->execute();

}