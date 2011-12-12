/** @file DetectorDisplay.cxx
@brief  Declare, implement class DetectorDisplay

* $Header$
*/

// includes
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/Property.h"

#include "GuiSvc/IGuiTool.h"
#include "gui/SubMenu.h"
#include "gui/GuiMgr.h"
#include "gui/Command.h"

#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"

#include "World.h"
#include "DisplayGeometry.h"
#include "MaterialDisplay.h"
#include <sstream>
namespace {
    gui::DisplayControl* pdisplay;
    const int nLevels = 9;
}
/** 
* @class DetectorDisplay
*
* @brief  A Tool that sets up the detector display 
* $Header$
*/
class DetectorDisplay : public AlgTool, virtual public IGuiTool 
{
public:

    //! Create the tool
    //! @param name The name of the service
    //!
    DetectorDisplay ( const std::string& type, const std::string& name, const IInterface* parent);

    virtual ~DetectorDisplay (){};

    /// implement to define gui elements: will be called from GuiSvc
    StatusCode initialize(gui::GuiMgr*);

private:
    // following specifically for the detector display
    FloatProperty m_axisSize;   //! axis size to use (in mm) 0 for no axis display
    IntegerProperty m_showLevel;    //! number of show commands to reveal details
    DisplayGeometry * m_glast;
#if defined(__GNUC__) && (__GNUC__ == 2)
    //---------- command class to set detail level-----------------
    class ShowDetector : public gui::Command {
    public:
        ShowDetector(gui::DisplayRep* rep, int level): m_rep(rep),m_level(level){}
        void execute(){ 
            m_rep->hide(false);    
            for(int i =0; i< m_level; ++i) m_rep->show(false); 
            pdisplay->redisplay();
        }
    private:
        gui::DisplayRep* m_rep;
        int m_level;
    };
    //------------------------------------------------
#endif
};


// declare the service factories for the DetectorDisplay
//static ToolFactory<DetectorDisplay> a_factory;
//const IToolFactory& DetectorDisplayFactory = a_factory;
DECLARE_TOOL_FACTORY(DetectorDisplay);

/// Standard Constructor
DetectorDisplay::DetectorDisplay(const std::string& type, 
                                 const std::string& name, 
                                 const IInterface* parent)
                                 : AlgTool( type, name, parent ) 
{

    // declare the properties
    declareProperty("axisSize", m_axisSize = 1000);
    declareProperty("showLevel", m_showLevel = 2);

    // Declare additional interface
    declareInterface<IGuiTool>(this);
}

// initialize
StatusCode DetectorDisplay::initialize (gui::GuiMgr* guiMgr) 
{
    StatusCode  status = StatusCode::SUCCESS;

    // bind all of the properties for this service
    status = setProperties ();

    // open the message log
    MsgStream log( msgSvc(), name() );

    // Get the Glast detector service, needed for constants, to keep the GlastDetector pointer
    IGlastDetSvc* gsv;
    if( service( "GlastDetSvc", gsv).isFailure() ) {
        log << MSG::ERROR << "Couldn't set up GlastDetSvc!" << endreq;
        return StatusCode::FAILURE;
    }

    // create the world instance, will be the top level of the tree 
    World* world = World::instance();

    log << MSG::DEBUG << "Instantiating the geometry using the GlastDetSvc " << endreq;
    m_glast = new DisplayGeometry();
    gsv->accept(*m_glast);
    if( log.level() <= MSG::DEBUG ) {
        log << MSG::DEBUG; 
        m_glast->printStats(log.stream());
        log << endreq;
    }
    // use the display from the Gui service to add the reps
    pdisplay = &guiMgr->display(); // put pointer into anonymous namespace
    pdisplay->setAxisSize(m_axisSize);

    // submenu for the dector
    gui::SubMenu& detmenu= guiMgr->display().menu().subMenu("Detector");
    gui::DisplayRep * worldRep= world->detectorViewer();
    worldRep->update();  // necessary to finish building nested reps

    // test for topLevel=LAT and presence of globalOffset_dz
    // this means that there's an extra level in the geometry, so
    // we bump the level counter  

    double temp;
    StatusCode sc = gsv->getNumericConstByName("globalOffset_dz", &temp);
    int thisLevel = 0;
    if (sc.isSuccess() && gsv->getTopVolume()=="LAT") thisLevel++;
    thisLevel += m_showLevel;
    if (thisLevel>nLevels) thisLevel = nLevels;

    for(int i =0; i< thisLevel; ++i) worldRep->show(false); 

    pdisplay->useMenu(&detmenu);
    pdisplay->add(worldRep, "detector", -1);
#if !defined(__GNUC__) || (__GNUC__ != 2)
    //---------- command class to set detail level-----------------
    class ShowDetector : public gui::Command {
    public:
        ShowDetector(gui::DisplayRep* rep, int level): m_rep(rep),m_level(level){}
        void execute(){ 
            m_rep->hide(false);    
            for(int i =0; i< m_level; ++i) m_rep->show(false); 
            pdisplay->redisplay();
        }
    private:
        gui::DisplayRep* m_rep;
        int m_level;
    };
    //------------------------------------------------
#endif

    for( int level=0; level<nLevels; ++level){
        std::stringstream label; label << "detail level " << level;
        detmenu.addButton( label.str(),  new ShowDetector(worldRep, level)); 
    }

    detmenu.addSeparator();
    // special display of all volumes with a given material
    MaterialDisplay mdxxx( guiMgr , m_glast, detmenu);
    pdisplay->useMenu();

    return status;
}



