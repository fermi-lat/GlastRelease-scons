
// Gaudi specific include files
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"

// gui, display includes
#include "GuiSvc/IGuiSvc.h"
#include "gui/DisplayControl.h"
#include "gui/GuiMgr.h"

/** @class AcdDisplay
 * @brief Display the ACD DOCA's as numbers at the tile positions
 *
 * @author Heather Kelly
 * $Header$
 */
class AcdDisplay : public Algorithm
{
public:
    // Constructor of this form must be provided
    AcdDisplay(const std::string& name, ISvcLocator* pSvcLocator); 
    virtual ~AcdDisplay() {}

    StatusCode initialize();
    StatusCode execute();
    StatusCode finalize(){ return StatusCode::SUCCESS;}
    
};

static const AlgFactory<AcdDisplay>  Factory;
const IAlgFactory& AcdDisplayFactory = Factory;

AcdDisplay::AcdDisplay(const std::string& name, ISvcLocator* pSvcLocator) :
Algorithm(name, pSvcLocator) 
{
    
}

/** @class AcdRep
 * @brief Handles the mechanics of updating the display
 *
 * @author Heather Kelly
 */
class AcdRep : public gui::DisplayRep {
public:
    AcdRep(){}
    void update(){
	// loop thru the tiles with DOCA's and show marker and value
#if 0 // example code from old version
	if( m_DOCA>100. ) return;
	set_color("red");
	markerAt(m_hit_tile.position());
	move_to(m_hit_tile.position()); 
	char text[10]; sprintf(text, "%5g", m_DOCA);
	drawText(text);
	set_color("black");
#endif // of example code 
    }
};

StatusCode AcdDisplay::initialize()
{
    //Look for the gui service
    IGuiSvc* guiSvc = 0;
    StatusCode sc = service("GuiSvc", guiSvc);

    
    //Ok, see if we can set up the display
    if (sc.isSuccess())  {
	//Set up the display rep for Clusters
	guiSvc->guiMgr()->display().add(new AcdRep(), "ACD recon");
    }
    
    return sc;
}

StatusCode AcdDisplay::execute()
{
    StatusCode sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());

    // set pointer to ACD data on TDS
    //    m_SiClusters = SmartDataPtr<SiClusters>(eventSvc(),"/Event/AcdRecon/SiClusters");
    
    return sc;
}


