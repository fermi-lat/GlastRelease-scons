/** 
* @file GuiSvc.h
* @brief declaration of the interface for the class GuiSvc, implementing IGuiSvc
*
*  $Header$
*/

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

/** @class GuiSvc
 * @brief Define a service that gives access to the GUI, specificially to a GuiMgr instance
 *
 * @author Toby Burnett
 * $Header$
 *
 * Implements Gaudi's IRunable interface in order to be able to control the event loop; for this requires
 * the following line in the job options file:
 *
 *   ApplicationMgr.Runable= "GuiSvc";
 *
 *  @section properties properties
 *
 *  The following properties are used to initialize the GuiMgr instance.
 *
 *     - size  [-300] Set initial size of the display window. Negative means z-axis upward
 *     - pause_interval [0] Time to wait, while displaying current event, before starting next
 *     - paused [true]  Initial setting
 *     - EvtMax [0xFEEDBABE] Internal use. Normally gets this property from the Application manager
 *
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
