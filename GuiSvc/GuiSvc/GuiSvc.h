// $Header$
// 
//  Original author: Toby Burnett tburnett@u.washington.edu

#ifndef _H_GuiSvc_
#define _H_GuiSvc_

// includes
#include "GaudiKernel/Service.h"
#include "GaudiKernel/IIncidentListener.h"
#include "GaudiKernel/IRunable.h"

#include "GuiSvc/IGuiSvc.h"

//forward declarations
namespace gui { class GuiMgr; }
template <class TYPE> class SvcFactory;
class IAppMgrUI;

/*!  
*/
class GuiSvc :  public Service, virtual public IIncidentListener, 
                virtual public IGuiSvc, virtual public IRunable
{  


public:
    /// only thing this wrapper class does: provide access to the wrapped class
    gui::GuiMgr* guiMgr() {return m_guiMgr;}

    void pause();

    // these needed to implement gui::MenuClient interface in order to schedule a quit call.
    void quit();
    void finishSetup(){};  // dummy

    //------------------------------------------------------------------
    //  stuff required by a Service
    
    /// perform initializations for this service. 
    virtual StatusCode initialize ();
    
    /// 
    virtual StatusCode finalize ();
    //------------------------------------------------------------------
    /// Handles incidents, implementing IIncidentListener interface
    virtual void handle(const Incident& inc);    
 
   /// Query interface
   virtual StatusCode queryInterface( const IID& riid, void** ppvUnknown );

   // for the IRunnable interfce
   virtual StatusCode run();


protected: 

    /// Standard Constructor
    GuiSvc ( const std::string& name, ISvcLocator* al );
    

private:
    // Allow SvcFactory to instantiate the service.
    friend class SvcFactory<GuiSvc>;
    void beginEvent();
    void endEvent();
    //
    void queryEvtMax();
    void queryPause();
    void queryOutputLevel();



    /// Reference to application manager UI
    IAppMgrUI*    m_appMgrUI;
    int m_evtMax;

    gui::GuiMgr* m_guiMgr;
    float   m_size;
    int     m_pause_interval;
    bool    m_paused;

};


#endif // _H_GuiSvc
