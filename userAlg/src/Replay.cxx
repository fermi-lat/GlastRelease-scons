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
  ApplicationMgr.TopAlg += {"Replay"};
  Replay.Members = {"TkrReconAlg", "TkrFindAlg"};
 </pre>

*/
class Replay : public Algorithm {
public:
    Replay(const std::string& name, ISvcLocator* pSvcLocator);

    StatusCode initialize();
    StatusCode execute();
    StatusCode finalize(){return StatusCode::SUCCESS;}
    
private:
    StatusCode addButton(const std::string & name);

    // Data Members

    StringProperty m_algname; 
    StringArrayProperty m_algnameList;
    IAlgManager* m_AlgMgr;

    // allow more than one instance to add to the same menu
    static gui::SubMenu* s_replayMenu;

};

static const AlgFactory<Replay>  Factory;
const IAlgFactory& ReplayFactory = Factory;

gui::SubMenu* Replay::s_replayMenu =0;

Replay::Replay(const std::string& name, ISvcLocator* pSvcLocator)
:Algorithm(name, pSvcLocator)
{
    // declare properties with setProperties calls
    declareProperty("algname", m_algname); // for backward compatibility
    declareProperty("Members",  m_algnameList);
    
}

//! set parameters and attach to various perhaps useful services.
StatusCode Replay::initialize(){
    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "initialize" << endreq;
    
    // Use the Job options service to set the Algorithm's parameters
    setProperties();
    
       
   // Alg manager service for finding algs by name
    if ( (serviceLocator( )->getService( "ApplicationMgr",
        IID_IAlgManager,
        (IInterface*&)m_AlgMgr )).isFailure() ) {
        log << MSG::ERROR << " no application manager???" << endreq;
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


    // add the button(s)
    if( m_algname.value() != "") sc = addButton(m_algname);

   std::vector<std::string> names = m_algnameList.value();
   for( std::vector<std::string>::const_iterator it=names.begin(); it!=names.end(); ++it){
       if( addButton(*it).isFailure() ) return StatusCode::FAILURE;
   }

   return sc;
}

StatusCode Replay::addButton(const std::string & algname)
{
    MsgStream log(msgSvc(), name());

    IAlgorithm * alg=0; 

    if( (m_AlgMgr->getAlgorithm(algname,  alg)).isFailure() ) {
        log << MSG::ERROR 
            << "Did not find the Algorithm "<< algname  << endreq;
        return StatusCode::FAILURE;
    }

    // define the recall class right here: all it does is call the alg.
    class ReplayCommand: public gui::Command {
    public:
        ReplayCommand(IAlgorithm* alg):m_alg(alg){}
        void execute(){
            m_alg->resetExecuted();
            m_alg->execute();
        }
    private:
        IAlgorithm * m_alg;
    };

    s_replayMenu->addButton(algname, new ReplayCommand(alg));
    return StatusCode::SUCCESS;

}

StatusCode Replay::execute()
{
    // execute method not really used or needed.
    return StatusCode::SUCCESS; 

}