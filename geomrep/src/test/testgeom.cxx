//  Simple main to test some features of the GUI interface

#include <string> // make sure gets loaded first
#include "gui/CompoundCommand.h"
#include "gui/GUI.h"
#include "gui/SimpleCommand.h"
#include "gui/PrintControl.h"
#include "gui/GuiMgr.h"

#include "geometry/Box.h"
#include "geometry/Cone.h"
#include "geometry/Ray.h"

#include "gui/DisplayControl.h"
#include "geomrep/ConeRep.h"

# include <strstream>

gui::DisplayControl* g_display;

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class Lines : public gui::DisplayRep {
public:
    
    
    Lines( const gui::GraphicsVector& p, const Volume& v )
        : m_p(p), m_v(v){
        using namespace gui;	
        // make a key update the display
        CompoundCommand* cmd = 
            new CompoundCommand(new SimpleCommand<DisplayControl>(g_display, &DisplayControl::clear));
        cmd->append(new SimpleCommand<DisplayControl>(g_display,&DisplayControl::update));
        Menu::instance()->register_key('c',cmd);
    }
    
    void draw_a_line(const Vector& n){
        Point p( m_p.x(),m_p.y(), m_p.z() );
        
        if( !m_v.inside(p) ){
            // starting outside
            setColor("green");
            move_to(p);
            double d = m_v.distanceToEnter(Ray(p, n), 100);
            if (d == FLT_MAX) {
                // missed, apparently
                d = 100;
                line_to(p+d*n);
                setColor("black");
                return;
            }else {
                // make a line to instersection
                Point p2 = p+d*n;
                line_to(p2);
                p = p2;
            }
        }
        // starting, or now inside - move to exit surface
        setColor("red");
        move_to(p);  // if not already
        double d = m_v.distanceToLeave(Ray(p, n), 100);
        Point p2 = p+d*n;
        line_to(p2);
        p = p2;
        
        // now outside -- finish up
        setColor("green");
        move_to(p);
        line_to(p+100*n);
        setColor("black");
    }
    
    
    void update() {
        Vector n(0,0,1);
        int count = 30;
        double delta = 2*M_PI/count;
        
        for(int i =0; i<count; ++i, n.rotateY(delta)){
            draw_a_line(n);
        }
        draw_a_line(Vector(0,1,0));
        draw_a_line(Vector(0,-1,0));
        
    }
    
private:
    const gui::GraphicsVector& m_p;
    const Volume& m_v;
};


//----------------------------------------------------------------------------
//              main

int main(int argc, char* argv[])
{
    
    // get pointer to the GuiMgr instance 
    gui::GuiMgr* guiMgr=gui::GuiMgr::instance();
    
    
    // get the display for below
    gui::DisplayControl & display = guiMgr->display();
    g_display = & display;
    
    
    // create object to test
    //Box box(50,50,50); 
    Cone cone(50, 5, 10, 30, 40);
    
    const Volume& vol =  cone;
    
    display.add(new ConeRep(cone), "object to test" );
    display.add(new Lines(display.reference_point(), vol), "lines to object");
    
    // start the message loop
    guiMgr->menu().run();
    
    return 0;
}
//===============================================================
void WARNING(const char* s){
    gui::GUI::instance()->inform(s);}
void FATAL(const char* s){
    gui::GUI::instance()->inform(s);}//cerr << s << endl; exit(-1);}
//===============================================================



