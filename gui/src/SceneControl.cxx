//     $Id$
//  Author: Toby Burnett
//

#include <string>
#include "gui/SceneControl.h"
#include "gui/ViewPort.h"
#include "gui/Scene.h"


#include <ctype.h>      // for tolower

namespace gui {

///////////////////////////////////////////////////////////////////////////////
// euler angles and perspective constant used as defaults for ViewPort objects
static float eulers[4][4] =
  { { 0.,  0.,  0., 0.},
    { 0., 90.,  0., 0.},
    {90., 90.,  0., 0.},
    {11.,110.,  4., 0.25} };
// set of euler angles with the z-axis up
static float eulers_z_up[4][4] =
  { { 0.,  0.,    0., 0.},
    {180., 90.,  90., 0.},
    {-90., 90.,  90., 0.},
    {-20.,-81., -90., 0.25} };


static const char  *title[]={
                "View 1, Front(X-Y)",
                "View 2, Side (Z-Y)",
                "View 3, Plan (X-Z)",
                "View 4, General" };

static const int  initialView = 1;             // Define initial selected view
static const int  initialSingle = 1;           // Define initial format

static GraphicsVector center2d(0.5,0.5,0);		// center for 2-d view

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SceneControl::SceneControl(ViewPort* projector, float scale, int initial_view )
//----------------------------------------------------------
// implement constructor for abstract base class to manage display
:  m_detail(0.01f)   // stuff for doCommand
,  m_axis(1)
{
    m_max_view = sizeof(title)/sizeof(char*);


    for( int i=0; i< m_max_view; i++) {
	if( scale>0)
	    m_views.push_back( &projector->copy(&m_reference_point, scale, eulers[i], title[i]));
	else
	    m_views.push_back( &projector->copy(&m_reference_point, -scale, eulers_z_up[i], title[i]));

    }

    select_view_number(scale>0?initial_view:2);
    set_single(initialSingle);

    //create a scene to control
    m_scene = new Scene;

    // setup 2-d scene and view
    m_2d_scene = new Scene;
    m_2d_view = &projector->copy(&center2d, 0.5, eulers[0], "" );

}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SceneControl::~SceneControl()
{
    delete m_scene;
    for( int i=0; i<m_max_view; i++)delete  m_views[i];
    delete m_2d_scene;
    delete m_2d_view;	
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void SceneControl::draw()
//-----------------------
{
    m_2d_view->set_quad(0,0);

    if( (m_scene->enabled()) ) { 
        
        if( m_single_view ){
            m_current_view->set_quad(0,0);
            m_scene->draw(m_current_view);
        }
        else{
            for( int i=0; i<m_max_view; i++){
                m_views[i]->set_quad(i+1, i==m_selected_view);
                m_scene->draw(m_views[i]);
            }
        }
    }
    m_2d_scene->draw(m_2d_view);

}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void SceneControl::redisplay(Draw2D* temp)
// if the ViewPorts are Projectors, redisplay using temporary canvas
// (i.e., to postcript)
{
    Draw2D* current = view().canvas();
    int i;
    for( i=0; i<m_max_view; i++) { m_views[i]->set_canvas(temp);}
    m_2d_view->set_canvas(temp);
    
    draw();

    for( i=0; i<m_max_view; i++) { m_views[i]->set_canvas(current);}
    m_2d_view->set_canvas(current);
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void SceneControl::resize(int width, int height)
{
    for(int i=0; i<m_max_view; i++)
	m_views[i]->setAspectRatio( (float)width/(float)height);
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void SceneControl::displaycurrent(int hl)
//--------------------------------------
{
    if(	m_single_view) m_current_view->set_quad(0,0);
    else m_current_view->set_quad(m_selected_view+1, hl);
    m_scene->draw(m_current_view);
    m_2d_scene->draw(m_2d_view);
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void SceneControl::select_view_number(int number)
//-----------------------------------------------
{
    if( number<1 || number>m_max_view ) return;
    m_current_view = m_views[number-1];
    m_selected_view = number-1;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void SceneControl::set_single(int value)
//--------------------------------------
{
    m_single_view= value;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void SceneControl::move_by(float x, float y)
//-----------------------------------------
//  adjust m_reference_point,in direction perpendicular to current view
{   

    // prevent motion if click out of view
    if( abs(x) < 1.0 && abs(y) < 1.0 )

	m_reference_point += m_current_view->transformToWorld(x,y);
} 
void SceneControl::move_left( float howfar){ move_by( howfar,0);}
void SceneControl::move_right(float howfar){ move_by(-howfar,0);}
void SceneControl::move_up(   float howfar){ move_by(0,-howfar);}
void SceneControl::move_down( float howfar){ move_by(0, howfar);}

//------------------------------------------------------------------
//                 adjust magnification (zoom)

void SceneControl::zoom_in(float howfar) { m_current_view->zoom(1+howfar);}
void SceneControl::zoom_out(float howfar){ m_current_view->zoom(1-howfar);}
void SceneControl::set_zoom(float z)     { m_current_view->setZoom(z);}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void SceneControl::reset_view()
//----------------------------
{   m_reference_point = m_origin; 
    m_current_view->reset();
}

//----------------------------------------------------------------
//             rotate viewer position

void
SceneControl::rotate_euler(const float* delta)
{
    m_current_view->rotateEuler(delta);
}

void
SceneControl::rotate_euler(float dx, float dy, float dz){
    float theta[3]={dx,dy,dz};
    m_current_view->rotateEuler(theta);
}

void SceneControl::rotate_x(float theta){ m_current_view->rotateX(theta);}
void SceneControl::rotate_y(float theta){ m_current_view->rotateY(theta);}
void SceneControl::rotate_z(float theta){ m_current_view->rotateZ(theta);}

//----------------------------------------------------------------
//             change perspective
void SceneControl::move_out(float d){   m_current_view->moveOut(d);}
void SceneControl::move_in(float d) {   m_current_view->moveIn(d);}
void SceneControl::set_view_distance(float d){ m_current_view->setPerspective(d);}

//----------------------------------------------------------------
//              set detail level
void SceneControl::set_detail(float d){ m_current_view->setDetail(d);}

void SceneControl::setViewPort (float xmin, float ymin, float xmax, float ymax)
{
  if (xmin > xmax) {   float h = xmax; xmax = xmin; xmin = h; }
  if (ymin > ymax) {   float h = ymax; ymax = ymin; ymin = h; }

   move_by ((xmin+xmax)/2, (ymin+ymax)/2);

   // zoom too much ?
   if (xmax-xmin <= 0.02  ||  ymax-ymin <= 0.02)
      return;

   // no, so zoom in
   if (xmax - xmin > ymax-ymin)
      zoom_in (2/(xmax-xmin) - 1);
   else
      zoom_in (2/(ymax-ymin) - 1);
}

//
///////////////////////////////////////////////////////////////////////////////
//                    Built-in character control stuff
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
int SceneControl::do_command(const char * cmd)
//--------------------------------------------
// doCommand -- for interactive modification of view parameters
//  return -1: unrecognized
//          0: ok, no redisplay needed
//          1: must redisplay quadrant -- displaycurrent(1);
//          2: redisplay all           -- display();
{
   int command = tolower(cmd[0]);
   int n = -1;

   float angle=0;
   switch( command ){
     case 0:
          // extended command
          command = cmd[1];
          switch( command ) {
            case 75: move_left(); return 1;
            case 72: move_up();   return 1;
            case 77: move_right(); return 1;
            case 80: move_down();  return 1;
            case 73: zoom_in();    return 1;
            case 81: zoom_out();   return 1;
            case 71: reset_view(); m_detail=0.01f; return 1;
            default: return -1;
          }
     case 'b': move_left();  return 1;        // Back
     case 'p': move_up();    return 1;         // uP
     case 'f': move_right(); return 1;         // Forward
     case 'n': move_down();  return 1;         // dowN
     case 'm': zoom_in();    return 1;         // Small
     case 's': zoom_out();   return 1;         // Large
     case 'r': reset_view();  m_detail=0.01f; return 1;
     case 'u': return 1;                   //Update

     case '0': n=0; break;
     case '1': n=1; break;
     case '2': n=2; break;
     case '3': n=3; break;
     case '4': n=4; break;

     case 'x': m_axis=0; return 0;
     case 'y': m_axis=1; return 0;
     case 'z': m_axis=2; return 0;
     case 'd': m_axis=3; return 0;

     case 'i': move_in(); return 1;
     case 'o': move_out(); return 1;

     case '>': angle=1.; break;
     case '.': angle=10; break;
     case '<': angle=-1; break;
     case ',': angle=-10; break;

     default: return -1;
   }
   if (angle ) {
                switch( m_axis ) {
              case 0: rotate_x(angle); return 1;
              case 1: rotate_y(angle); return 1;
              case 2: rotate_z(angle); return 1;
              case 3: set_detail(m_detail *= angle>0? 2. : 0.5 ); return 1;
           }
   }

    // here if view changed

    if( n==m_selected_view+1 || n==-1) return 0;
    if (n==0) {
       // toggling quad views 
         if (m_single_view)  set_single(0);
         else            set_single(1);
         return 2;
    }
    select_view_number(n);
    return 2;

}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void 
SceneControl::print_commands(std::ostream& cout)
{
cout <<
"\n"
"Keystroke display commands:\n"
"===========================\n"
"  b,f  \t move left,right  (Back,Forward)\n"
"  p,n  \t move up,down     (uP, Down)\n"
"  m,s  \t zoom in,out      (Magnify, Shrink)\n"
"  i,o  \t move in,out      (In, Out)\n"
"  r    \t reset view  (Reset)\n"
"  0    \t toggle single/quad view\n"
"  1,2,3,4\t select  view#\n"
"  u    \t update     (Update)\n"
"\n"
"Following set action for ,.\n"
"  x,y,z\t set x,y,z axis for rotation\n"
"  d    \t modify detail level by x2 \n"
"  ,    \t (comma)  rotate left by 10 deg\n"
"  <    \t rotate left by 1 deg\n"
"  .    \t (period) rotate right by 10 deg\n"
"  >    \t rotate right by 1 deg\n"
;
    cout.flush();
}

} // namespace gui


