//     $Id$
//  Author: Toby Burnett
//
//

#ifndef DISPLAYCONTROL_H
#define DISPLAYCONTROL_H


#include <string>
#include <vector>

#include "gui/Command.h"
#include "gui/Menu.h"

#include "gui/DisplayRep.h"
#include "gui/DisplayRep2D.h"

namespace gui {

class Scene;
class SceneControl;
class Draw2D;

class DisplayControl  {
    // definition of the class DisplayControl--handle event display. It
    // has a list of <<DisplayRep>> objects which it manages, and sets up a "Display"
    // menu.

public:

    DisplayControl(Menu& menu, SceneControl* control); 
    // constructor with reference to menu, scene controler

    ~DisplayControl();
    
    GUI::Toggle* add(DisplayRep* arep, const std::string& name=std::string(""), int state=1);
    // Add an element to the display
    // It the name argument is present, create a  button for it under the "display"
    // pull-down menu. If state is 0 or 1 it will be a toggle, starting show/hide. 
    // In that case, return pointer to a Toggle that allows independent set/unset
    // Otherwise create separate show and hide buttons

    void add(DisplayRep2D* arep, const std::string& name=std::string(""), int state=1);
    // add  element to 2-d display

    void setTitle(const std::string& title);
    // set the title bar of the associated window

    void setAxisSize(float size=1);

    GraphicsVector& reference_point();
    const GraphicsVector& reference_point()const;
    // access to current reference point
    
    void postScript();
    // dump current view to a file for postscript (if viewport is a Projector)

    void vrml(std::ostream& out);
    // write out a VRML file

    void update();
    // ask scene to update

    void redisplay();
    // redisplay all objects, if running

    void set_running(bool s){m_running=s;}
    // use to temporarily disable redisplay when modifying show/hide bit of many objects
	 

    void clear();
    // clear display (send message to the scene)
    
    Scene& scene();
    const Scene& scene()const;
    // access to the 3-d scene

    SubMenu& menu(); 
    // return the (sub)menu associated with the display

    void useMenu(SubMenu* s=0);
    // set a different submenu (which could be a sub-submenu of this Display menu)
    // if no argument, revert to the standard

    // needed to implement interface MenuClient
    void finishSetup();    // adds menu items after all user buttons
    void quit(){};

    //------------------------------------------------------------------
    // Define nested class to handle submenu
    class SubMenu     {
        // This is a special SubMenu, containing a pointer to a gui::SubMenu,
        // and a list of nested SubMenus (allowing a tree)
    public:
        ~SubMenu();
        
        SubMenu& subMenu(const std::string& label, DisplayRep* rep=0);
        // return reference to a new sub menu, with given button label, and optional rep
        // that will be controlled by a "Show All" button at its top.
        // This button also applies to all submenus of this

        void add(DisplayRep * rep, const std::string& name="");
        // add a rep with optional name to be used for button

    private:
        friend class DisplayControl;
        void hide();
        void hide(bool update);
        void show();
        void show(bool update);
        SubMenu(DisplayControl* display, SubMenu* parent, DisplayRep * rep, const std::string& name);
        DisplayRep* _rep;  // optional rep
        gui::SubMenu* _menu; // 
        DisplayControl* _display;
        std::vector<SubMenu*> _submenus;
        std::vector<GUI::Toggle*> _rep_list; // list of toggles for reps
    };
    //------------------------------------------------------------------

    SubMenu& subMenu(const std::string& label, DisplayRep* rep=0);
    // function that returns a new SubMenu


private:
    void setRefPt(); // called by registered button to get dialog
    void setScale(); // called by registered button to get dialog
    GUI::Toggle* addButton(DisplayRep* arep, const std::string& name, int state);
    // internal routine
    
    gui::SubMenu*  m_sub_menu; // special sub menu we control
    gui::SubMenu*  m_user_menu;// user-settable menu
    Menu*	m_menu_bar;     // original menu object

    Scene*	m_scene;   	// Scene (3d )
    Scene*      m_2d_scene;     // 2d scene
    SceneControl* m_control;	// Scene controller
    
    float 	m_axis_size;
    bool	m_running;	// used to keep track of startup

    Draw2D*     m_canvas;	// access to the 2-d canvas

    std::vector<SubMenu*> _submenus; // list to delete.
};

inline gui::SubMenu& DisplayControl::menu(){return *m_sub_menu;} 
inline gui::Scene&   DisplayControl::scene(){return *m_scene;}
inline const gui::Scene&   DisplayControl::scene()const{return *m_scene;}

} // namespace gui
#endif 


