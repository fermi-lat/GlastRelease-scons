//
// Gaudi system includes
#include "GaudiKernel/Algorithm.h"

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IAlgManager.h"

// gui includes
#include "GuiSvc/IGuiSvc.h"
#include "gui/SubMenu.h"
#include "gui/GuiMgr.h"
#include "gui/Command.h"
#include "gui/Menu.h"


/** @class Replay 
 * @brief Call (replay?) a given algorithm by name, settable from job options
 * 
 * $Header$
 * <br> Following is example of how to select multiple algs:
<pre>
  ApplicationMgr.TopAlg += {"Replay/button1", "Replay/button2"};
  button1.algname="TkrReconAlg";
  button2.algname="TkrFindAlg";
 </pre>

*/
class Replay : public Algorithm {
public:
    Replay(const std::string& name, ISvcLocator* pSvcLocator);

    StatusCode initialize();
    StatusCode execute();
    StatusCode finalize(){return StatusCode::SUCCESS;}
    
private:
    // Data Members

    IAlgorithm* m_alg;
    StringProperty m_algname;
    static gui::SubMenu *s_replayMenu;

};

static const AlgFactory<Replay>  Factory;
const IAlgFactory& ReplayFactory = Factory;

gui::SubMenu* Replay::s_replayMenu =0;

Replay::Replay(const std::string& name, ISvcLocator* pSvcLocator)
:Algorithm(name, pSvcLocator)
{
    // declare properties with setProperties calls
    declareProperty("algName",  m_algname);
    
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

    // get the  Gui service, required for interactive function
    //
    IGuiSvc* guiSvc=0;
    
    if ( service("GuiSvc", guiSvc).isFailure() ){
        log << MSG::WARNING << "No GuiSvc, no display" << endreq;
        return StatusCode::SUCCESS;  // just ignore thisy     
    }
    
    // found: create new top-level menu
    gui::GuiMgr* guiMgr=guiSvc->guiMgr() ;
    if( s_replayMenu==0) {
        s_replayMenu = & guiMgr->menu().subMenu("Replay");
    }

    // define the recall class right here: all it does is call the alg.
    class ReplayCommand: public gui::Command {
    public:
        ReplayCommand(IAlgorithm* alg):m_alg(alg){}
        void execute(){
            m_alg->execute();
        }
    private:
        IAlgorithm * m_alg;
    };
    s_replayMenu->addButton(m_algname, new ReplayCommand(m_alg));

   return sc;
}

StatusCode Replay::execute()
{
    // execute method not really used or needed.
    return StatusCode::SUCCESS; 

}