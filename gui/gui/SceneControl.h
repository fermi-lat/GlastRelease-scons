//     $Id$
// Project: Arve graphics
//  Author: Toby Burnett
//

#ifndef  SCENECONTROL_H
#define  SCENECONTROL_H

#include <iostream>

#include "gui/ViewPort.h"
#include <vector>

namespace gui {

class Scene;
class Draw2D;

class SubMenu;


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class SceneControl 
//  *Manage an Scene.    
//  *Has a list of (now four) ViewPort objects, which know how to
//    draw contents of DisplayRep's 
//  *Provides interface to ViewPort projection options
//  *Implements keystroke commands to modify views
// Also manage a second Scene, with 2-d transform 
{
public:

    
    virtual ~SceneControl();
    
    GraphicsVector & reference_point(){return m_reference_point;}
    // get a modifiable reference to the current display refrence point

    void draw();
    // generate a display

    virtual void redisplay(){draw();}
    // subclass that knows about display device can implement to
    // trigger refresh

    void redisplay(Draw2D* temp);
    // if the ViewPorts are Projectors, redisplay using temporary canvas
    // (i.e., to postcript)

    void displaycurrent(int hl);
    // redisplay the current window

    void select_view_number(int number);
    // select view: number=1-4 for view number

    void set_origin(const GraphicsVector& neworig){m_origin = neworig;}
    // modify basic display

    void set_single(int value);
    // select view and display it only (0: means all)

    void reset_view();
    // reset reference point, current view parameters

    ViewPort& view();
    // return reference to current view

    Scene& scene()const;
    // reference to the scene object

    Scene& scene_2d()const;

    //-----------------------------
    // adjustments to current view
    //-----------------------------

    void move_by(float,float);
    void move_left(float x=0.1);
    void move_right(float x=0.1);
    void move_up(float x=0.1);
    void move_down(float x=0.1);
    // adjust reference point according to current view

    void zoom_in(float x = 0.25);
    void zoom_out(float x= 0.25);
    void set_zoom(float z = 1.);
    // adjust zoom factor

    void move_in(float=0.1);
    void move_out(float=0.1);
    void set_view_distance(float=0);
    // adjust perspective

    void set_detail(float d=0.);
    // set detail threshold

    void setViewPort (float xmin, float ymin, float xmax, float ymax);
    // pan and zoom to the given boundaries

    void rotate_x(float);
    void rotate_y(float);
    void rotate_z(float);
    void rotate_euler(const float*);
    void rotate_euler(float,float,float);
    // change orientation of current view

    int do_command(const char *);
    // display control, allowing interactive modification of view parameters
    //  return -1: unrecognized
    //          0: ok, no redisplay needed
    //          1: must redisplay current view
    //          2: must redraw entire screen
    // Subclass can implement a displayer by passing keystroke commands to
    // it in an enent loop


    virtual void print_commands(std::ostream& );
    // prints list of commands to the stream

    void set_menu(const SubMenu* m){m_menu = m;}
    const SubMenu* menu()const{return m_menu;}

protected:
    SceneControl(ViewPort*, float size=100 , int initial_view=1);
    // constructor, for subclass that will fill in the ViewPort object

    void resize(int x, int y);	
    // process screen size change

    
private:

    Scene*	m_scene;
    Scene*	m_2d_scene;
    // the scenes that we control, 3d multiple views, and 2d

    GraphicsVector m_origin;                // default reference point
    GraphicsVector m_reference_point;        // current reference point

    int	    m_max_view;		    // number of views
    int	    m_single_view;          // Flag if single or multiple views
    int	    m_selected_view;	    // selected view number if multiple displayed


    std::vector<ViewPort*> m_views;	    // list of active viewports
    ViewPort*         m_current_view;       //  current view
    ViewPort*	    m_2d_view;		// special view for 2-d display


    float m_detail; // doCommand
    int   m_axis;
    bool	m_enabled;	// toggle for enabling display

    const SubMenu* m_menu;    // pointer to associated menu


};
//----------------------- inlines ---------------------------------
inline ViewPort&
SceneControl::view(){
    return *m_current_view;
}
inline Scene& SceneControl::scene()const{return *m_scene;}
inline Scene& SceneControl::scene_2d()const{return *m_2d_scene;}

} // namespace gui

#endif

