//     $Id$
//  Author: Toby Burnett
//
// implementation of  Display control

#ifdef NT_MSVCPP
#pragma warning(disable: 4786)
#endif

#include "gui/DisplayControl.h"
#include "gui/GUI.h"   // for the gui window
#include "gui/SubMenu.h"

#include "gui/DisplayRep.h"
#include "gui/SceneControl.h"
#include "gui/Scene.h"
#include "Vrml.h"
#include "gui/Draw2D.h"
#include "PSdraw.h"
#include "gui/Projector.h"


#include "gui/PrintControl.h"
#include "gui/SimpleCommand.h"

#include <string>  
#include <fstream>
#include <strstream>
#include <iomanip>
#include <cmath>

namespace gui {

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
DisplayControl::DisplayControl(Menu& menu, SceneControl* control)
    : m_control(control),
     m_axis_size(100),
     m_running(false)
{
  m_menu_bar = &menu;  // pointer to use below
    m_user_menu =m_sub_menu = &menu.subMenu("Display");

    // schedule to have our finishSetup called
    menu.add(new MenuClient<DisplayControl>(this));


    // get a pointer to the scene
    m_scene = &m_control->scene();
    m_2d_scene = &m_control->scene_2d();

    m_canvas = m_control->view().canvas();

    // give the scene control a pointer to the menu, so it can pop it up
    m_control->set_menu(m_sub_menu);


    // this is to force link to include the Projector
    Projector p((Draw2D*)0); 

}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
DisplayControl::~DisplayControl()
{
    delete m_control;
    for( std::vector<DisplaySubMenu*>::iterator it=_submenus.begin(); it !=_submenus.end(); ++it)
        delete *it;
}
//=====================================================================
//         anonymous namespace for helper classes
//=====================================================================
#ifdef NT_MSVCPP  //KAI can't handle this???
namespace {
#endif
double sqr(double x){return x*x;}

//     VRML_command class
class VRML_command : public Command {
public:
    VRML_command(DisplayControl* dcon):m_dcon(dcon){}
    void execute(); 
private:
    DisplayControl* m_dcon;
};

void VRML_command::execute() {
    static std::string filename("glast.wrl");
    char* answer = GUI::instance()->askForFileName(filename.c_str(),"save"
				,"vrml files (*.wrl)|*.wrl|all files|*.*|") ;
    if( !answer) return;	// check for cancel
    filename = answer;
    std::ofstream vrml_file(filename.c_str());
    if( vrml_file.rdbuf()->is_open() ) m_dcon->vrml(vrml_file);
    else GUI::instance()->inform("Could not open the file");
}

class PrintInstructions : public Command {
    friend class DisplayControl;
    PrintInstructions(SceneControl* disp, std::ostream* out):m_disp(disp), m_out(out){}
    void execute(){m_disp->print_commands(*m_out);}
    SceneControl* m_disp;
    std::ostream* m_out;

};

//   Private class used by DisplayControl to define axes
class Axes: public DisplayRep
{
public:
    Axes(double scale=1.);
    // the constructor creates the representation

    // dummy  clear, update
    void clear(){};
    void update(){};
    
private:
};


Axes::Axes(double scale)
{
    static GraphicsVector unitvector[3] =
    {GraphicsVector(1.,0.,0.),
    GraphicsVector(0.,1.,0.),
    GraphicsVector(0.,0.,1.)
    };
    
    static const char* axislabel[3]={"x","y","z"};
    set_color("black");
    for( int i=0; i<3; i++){
        move_to(GraphicsVector(0,0,0));
        line_to((scale)*unitvector[i]);
        drawText(axislabel[i]);
    }
    flush();
}

//=====================================================================
//   Private class used by DisplayControl to display reference point
static double  amin(double x){return  fabs(x)<1e-5? 0:x;}
class ReferencePoint : public DisplayRep2D
{
public:
    ReferencePoint(DisplayControl& dc):m_dc(dc){};
    void draw2D(Draw2D* canvas) {
        std::strstream text;
        GraphicsVector p(m_dc.reference_point());
        static int precision=3, w=10;
        text << "ref pt: "  << std::setprecision(precision) << std::setiosflags(std::ios::fixed)
	     << "("
	     << std::setw(w) << amin(p.x())  << ','
	     << std::setw(w) << amin(p.y())  << ','
	     << std::setw(w) << amin(p.z())  << ')' << '\0';
	const float s=0.025f;
	canvas->set_rgb(1,0,0);
	canvas->move_to(0,s);
	canvas->line_to(0,-s);
	canvas->move_to(s,0);
	canvas->line_to(-s,0);
	canvas->set_rgb(0,0,0);

	// display in ULH corner for now?
        canvas->draw_string(-0.95f, 0.9f, text.str());
#ifdef _MSC_VER // avoids memory leak in Visual
        text.freeze(false);
#endif
    }
private:
    const DisplayControl& m_dc;
};
//=====================================================================
//   Private class used by DisplayControl to display reference point
class Scale : public DisplayRep2D
{
public:
    Scale(SceneControl* dc):m_dc(dc){};
    void draw2D(Draw2D* canvas){
        const float xa=-0.9f, xb=-0.7f, y=0.8f;
        double scale =   
            sqrt(( m_dc->view().transformToWorld(xa,y)
                 - m_dc->view().transformToWorld(xb,y)
                 ).mag2());
        std::strstream text;
        text <<  std::setprecision(3) << scale  << '\0';
        canvas->move_to(xa, y+0.025);
        canvas->line_to(xa, y-0.025);
        canvas->move_to(xa, y);
	canvas->line_to(xb, y);
        canvas->move_to(xb, y+0.025);
	canvas->line_to(xb, y-0.025);
        canvas->draw_string(xb+0.01, y-0.01, text.str());

#ifdef _MSC_VER // avoids memory leak in Visual
        text.freeze(false);
#endif
    }
private:
    SceneControl* m_dc;
};
///////////////////////////////////////////////////////////////////////////////
//   add a view to display, with a buttons or a toggle in the Display menu
/////////////////////////////////////////////////////////////////////////////////
// File-scope stuff to define hide and show commands
    class ShowIt : public Command    {
    friend class DisplayControl;
	ShowIt( DisplayControl* s,DisplayRep * rep, bool f)
	    :m_scene(s),m_rep(rep),m_flag(f){}
	void execute()
	{
	    m_rep->show(m_flag);
	    m_scene->redisplay();
	}
	DisplayControl* m_scene;
	DisplayRep* m_rep;
	bool m_flag;
    };
    class HideIt : public Command    {
    friend class DisplayControl;
	HideIt( DisplayControl* s,DisplayRep * rep, bool f)
	    :m_scene(s),m_rep(rep),m_flag(f){}
	void execute()
	{
	    m_rep->hide(m_flag);
	    m_scene->redisplay();
	}
	DisplayControl* m_scene;
	DisplayRep* m_rep;
	bool m_flag;
    };

#ifdef NT_MSVCPP
} // anonymous namespace
#endif
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void
DisplayControl::finishSetup()
//---------------------------
{
    m_sub_menu->addSeparator();

    if( m_axis_size != 0)  add( new Axes(m_axis_size),"axes" );
    add( new ReferencePoint(*this),"Reference point", false);
    add( new Scale(m_control), "Scale", false);

    m_sub_menu->addButton("Set ref. pt...", 
        new SimpleCommand<DisplayControl>(this, &DisplayControl::setRefPt));

    m_sub_menu->addButton("Set magnification...", 

        new SimpleCommand<DisplayControl>(this, &DisplayControl::setScale));

    // the rest of these on a sub menu

    gui::SubMenu* system= &m_sub_menu->subMenu("System");

    system->addToggle("3D Enabled", m_scene->enabled(),
	new SimpleCommand<Scene>(m_scene,&Scene::enable),
	new SimpleCommand<Scene>(m_scene,&Scene::disable)
		);
    system->addToggle("2D Enabled", m_2d_scene->enabled(),
	new SimpleCommand<Scene>(m_2d_scene,&Scene::enable),
	new SimpleCommand<Scene>(m_2d_scene,&Scene::disable)
		);
    system->addButton("clear", 
	new SimpleCommand<DisplayControl>(this,&DisplayControl::clear));

    system->addButton("update", 
	new SimpleCommand<DisplayControl>(this,&DisplayControl::update));

    system->addButton("Dump to postscript",
	new SimpleCommand<DisplayControl>(this,&DisplayControl::postScript));
    system->addButton("Dump to vrml",new VRML_command(this));
    system->addButton("Print instructions",
			   new PrintInstructions(m_control,PrintControl::instance()->out()));

    m_running = true;
    m_scene->update();
    m_2d_scene->update();
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void DisplayControl::setRefPt()  // called by button
{
   GraphicsVector&  p = m_control->reference_point();  
   double d[]={p.x(),p.y(),p.z()};
   m_menu_bar->query("Reference point?",d,3);
   p = GraphicsVector(d[0],d[1],d[2]);
 
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void DisplayControl::setScale()  // called by dialog 
{
    double scale = m_control->view().zoom(); //current zoom factor
    double zoom_factor = scale;
    m_menu_bar->query("Magnification?",&scale);
    if(zoom_factor==scale&&zoom_factor!=0) return;
    m_control->view().zoom(scale/zoom_factor); 
    m_control->redisplay();
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void DisplayControl::setAxisSize(float size)
{
    m_axis_size = size;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void DisplayControl::redisplay()
//---------------------------
// This is necesary to prevent the setup calls for ShowIt and HideIt calling
// Motif before the display has been realized
{
    if( m_running )
	m_control->redisplay();
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void DisplayControl::clear()
//-------------------------
{
    m_scene->clear();
    m_2d_scene->clear();
    redisplay();
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void DisplayControl::update()
//-------------------------
{
    m_scene->update();
    m_2d_scene->update();
    redisplay();
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
GUI::Toggle* 
DisplayControl::add(DisplayRep* pRep, const std::string& name, int state)
{
    m_scene->add(pRep);
	
    // now control it, if requested
    if( !name.length() ) return 0;

    return addButton(pRep,name,state);
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void
DisplayControl::add(DisplayRep2D* pRep, const std::string& name, int state)
{
    m_control->scene_2d().add(pRep);

    if( name.empty() ) return;
    addButton(pRep, name, state);

}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
GUI::Toggle*
DisplayControl::addButton(DisplayRep* pRep, const std::string& name, int state)
{

    Command* showit = new ShowIt(this, pRep, state>=0);
    Command* hideit = new HideIt(this, pRep, state>=0);
    if( state <0 ) {
	// special flag to make separate show and hide buttons
	m_user_menu->addButton(std::string("show ")+name, showit);
	m_user_menu->addButton(std::string("hide ")+name, hideit);
        return 0;
    } 

    if( state>0 ) showit->execute(); else hideit->execute();
    return m_user_menu->addToggle(name, state>0, showit, hideit );

}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void DisplayControl::useMenu(gui::SubMenu* m)
{
    m_user_menu = m? m : m_sub_menu;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void
DisplayControl::setTitle(const std::string& t){GUI::instance()->setTitle(t.c_str());}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
GraphicsVector&  DisplayControl::reference_point()
//--------------------------------------------------
{
  return m_control->reference_point();
}

const GraphicsVector&  DisplayControl::reference_point()const
//--------------------------------------------------
{
  return const_cast<const GraphicsVector&>(m_control->reference_point());
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void DisplayControl::vrml(std::ostream& pout)
//--------------------------------------
// ask the scene to dump itself to a VRML file
{
  Vrml output(pout);
    ((Scene*)m_scene)->draw(&output);
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void DisplayControl::postScript()
//--------------------------------
// ask the scene to draw itself to a file in postscript
{
    static const char* filename = "temp.ps";
    PSdraw ps(filename);
    if( !ps() ) {
	GUI::instance()->inform("Cannot open file for postscript output");
	return;
    }

    m_control->redisplay(&ps);
    std::ostrstream msg; msg << "Postscript version written to file " << filename << '\0';
    GUI::instance()->inform(msg.str());
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//      SubMenu implementation
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
DisplayControl::DisplaySubMenu& DisplayControl::subMenu(const std::string& name, DisplayRep *rep)
{
    DisplaySubMenu* s = new DisplaySubMenu(this, 0, rep, name);
    _submenus.push_back(s);  // save in this list to delete 
    return *s; 
}

DisplayControl::DisplaySubMenu::DisplaySubMenu(DisplayControl* display,DisplaySubMenu* parent, DisplayRep * rep, const std::string& name)
: _rep(rep)
, _display(display)
{

    // create a new submenu with the name;
    if ( parent==0) {
        _menu = &_display->menu().subMenu(name) ;
    } else {
        _menu = &parent->_menu->subMenu(name); 
        // and add to list of DisplayReps in parent menu
        parent->_submenus.push_back(this);
    }

    // add the rep, if there is one, to the DisplayControl object for display
    if( rep) _display->add(rep);

    // top buttons to show or hide all on the menu or under it
    _menu->addButton("Show all", new gui::SimpleCommand<DisplaySubMenu>(this, &DisplaySubMenu::show));
    _menu->addButton("Hide all", new gui::SimpleCommand<DisplaySubMenu>(this, &DisplaySubMenu::hide));
    _menu->addSeparator();


}
DisplayControl::DisplaySubMenu& DisplayControl::DisplaySubMenu::subMenu(const std::string& name,DisplayRep * rep)
{
    return *new DisplayControl::DisplaySubMenu(_display,this,rep,name);
}
    
void DisplayControl::DisplaySubMenu::add(DisplayRep * rep, const std::string& name)
{
    _display->useMenu(_menu);
    GUI::Toggle* display_toggle = _display->add(rep,name);
    if( display_toggle !=0) _rep_list.push_back(display_toggle);
    _display->useMenu();
}
void DisplayControl::DisplaySubMenu::show(){
   _display->set_running(false); show(true); _display->set_running(true); _display->redisplay();}
void DisplayControl::DisplaySubMenu::hide(){
   _display->set_running(false); hide(true); _display->set_running(true); _display->redisplay();}

void DisplayControl::DisplaySubMenu::show(bool update)
{
    if(_rep) _rep->show();
    for(std::vector<GUI::Toggle*>::iterator itt=_rep_list.begin(); itt!=_rep_list.end(); ++itt)
        (*itt)->set(); //set the toggle to the correct state
    for(std::vector<DisplaySubMenu*>::iterator it = _submenus.begin(); it!=_submenus.end(); ++it)
        (*it)->show(false);
}
void DisplayControl::DisplaySubMenu::hide(bool update)
{
    if(_rep) _rep->hide();
    
    for(std::vector<GUI::Toggle*>::iterator itt=_rep_list.begin(); itt!=_rep_list.end(); ++itt)
        (*itt)->unset();
    for(std::vector<DisplaySubMenu*>::iterator it = _submenus.begin(); it!=_submenus.end(); ++it)
        (*it)->hide(false); 
}

DisplayControl::DisplaySubMenu::~DisplaySubMenu()
{
    for(std::vector<DisplaySubMenu*>::iterator it = _submenus.begin(); it!=_submenus.end(); ++it) delete (*it);

}




}// namespace gui



