
// Gaudi specific include files
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"

#include "Event/TopLevel/EventModel.h"
#include "Event/Recon/AcdRecon.h"

#include "geometry/Point.h"

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
    AcdRep(IDataProviderSvc* dps){}
    void update(){
        SmartDataPtr<Event::AcdRecon> acdRec(m_dps, EventModel::AcdRecon::Event);
        if (!acdRec) return;

	// loop thru the tiles with DOCA's and show marker and value
	if( acdRec->getDoca()>1000. ) return;
	//set_color("red");
	//markerAt(Point(0,0,0));
	//move_to(Point(0,0,0)); 
	//char text[10]; sprintf(text, "%5g", acdRec->getDoca());
	//drawText(text);
	//set_color("black");
    }
private:
    IDataProviderSvc* m_dps;
};

StatusCode AcdDisplay::initialize()
{
    //Look for the gui service
    IGuiSvc*   guiSvc = 0;
    StatusCode sc     = service("GuiSvc", guiSvc);
    if( sc.isFailure() )  return sc;

    
    //Ok, see if we can set up the display
    if (sc.isSuccess())  {
        gui::DisplayControl& display = guiSvc->guiMgr()->display();
        
        gui::DisplayControl::DisplaySubMenu& acdmenu = display.subMenu("AcdRecon");

	guiSvc->guiMgr()->display().add(new AcdRep(eventSvc()), "ACD recon");
    }
    
    return sc;
}

StatusCode AcdDisplay::execute()
{
    return StatusCode::SUCCESS;
}


