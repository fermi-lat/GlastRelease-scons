// MaterialDisplay.cxx: implementation of the MaterialDisplay class.
//
// $Header$
//  Author: T. Burnett
//////////////////////////////////////////////////////////////////////

#include "MaterialDisplay.h"
#include "CompositeMedium.h"
#include "World.h"
#include "DisplayGeometry.h"

#include "geometry/Box.h"
#include "geometry/Tube.h"
#include "geomrep/BoxRep.h"
#include "geomrep/TubeRep.h"

#include "gui/Command.h"
#include "gui/GuiMgr.h"


#include <typeinfo>

static std::string color("red"); // current color to use for display
static gui::Menu* theMenu=0;       // set for use below
//============================================================================
//    GUI interface class 
class MaterialDisplay::ChooseColor : public gui::Command {
friend class MaterialDisplay;
    void execute(){
        theMenu->query("Select a color for next display", &color);
    }
};
//============================================================================
class MaterialDisplay::Rep : public gui::DisplayRep {
public:
    Rep( const std::string& material):m_material(material){};
    virtual ~Rep(){}
    void update();
    void clear(){};
private:
    void doit(const Medium* med);
    void innerLoop(const Medium* med);
    const std::string & m_material;
};
//============================================================================
MaterialDisplay::MaterialDisplay(gui::GuiMgr* guiMgr, DisplayGeometry* geom,
                                 gui::SubMenu& detmenu )
                                 :m_display(&(guiMgr->display())) 
{
    theMenu = &(guiMgr->menu());

    // make a submenu in the Display|Detector   submenu, tell display to use it for now
    gui::SubMenu& m = detmenu.subMenu("Materials");
    m_display->useMenu( &m  ); 

    m.addButton("Set color ...", new ChooseColor);
    m.addSeparator();
    // loop over known Materials
    const DisplayGeometry::MaterialSummary& materials = geom->materials();
    for(  DisplayGeometry::MaterialSummary::const_iterator  imat = materials.begin(); imat!=materials.end(); ++imat){
        const std::string& name = imat->first;
        if( name.empty() || name=="vacuum") continue;
        m_display->add(new Rep( imat->first), imat->first, false); 
    }
    m_display->useMenu(); // restore previous display menu
}
MaterialDisplay::~MaterialDisplay()
{
    m_display->useMenu(); // make sure reset
}
//============================================================================
//    MaterialDisplay::Rep implementation 

void MaterialDisplay::Rep::doit(const Medium* med)
{
    //  get the associated volume, check its type
    const Shape& vol = med->volume();
    const std::type_info& t = typeid(vol);
    if( t==typeid(Box) )    append(BoxRep(vol ));
    else if( t==typeid(Tube) )    append(TubeRep(vol ));
    else {
        std::cerr << "no display representation was set for volume " << t.name()  << std::endl; 
    }
}

void MaterialDisplay::Rep::innerLoop(const Medium* med)
{
    if( med->material() == m_material) doit(med);  // process this medium
    if(!med->isComposite())return;

    // it is composite: recurse
    const CompositeMedium* cmed = dynamic_cast<const CompositeMedium*>(med);
    for( CompositeMedium::const_iterator it = cmed->begin(); it < cmed->end(); ++it){
        innerLoop(*it);
    }
}

void MaterialDisplay::Rep::update()
{
    if( hasRepresentation())return; // only do this once
    setColor(color); // use color that user may have set
    innerLoop(World::instance());
    setColor("black"); // restore to default
}


