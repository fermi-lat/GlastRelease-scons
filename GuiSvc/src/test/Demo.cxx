// 
#include "Demo.h"

// gui event graphics include files
#include "gui/GuiMgr.h"
#include "gui/DisplayControl.h"
#include <strstream>
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// A class with 3-d info. It will be represented by a cube beween (a,a,a) and (b,b,b),
// and has a name that should plotted in the middle. 
class SimpleBox 
{
public:
    SimpleBox(double a, double b, const std::string& name)
        :m_a(a), m_b(b), m_name(name){}
    ~SimpleBox(){}
    double a()const{return m_a;}
    double b()const{return m_b;}
    const std::string& name()const{return m_name;}
private:
    double m_a, m_b;
    std::string m_name;
};
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// The class that makes a representation for plotting. It must implement
// an update() method to create the representation
class SimpleBoxRep : public gui::DisplayRep
{
public:
    SimpleBoxRep(const SimpleBox* pbox):m_pbox(pbox)
    {
       // update(); // make rep in constructor
    }
    void update()
    {
        using namespace gui;
        if( hasRepresentation() ) return;
        set_color("blue");
        double a = m_pbox->a(), b = m_pbox->b();

        // create a list of the vertices of the box
        addVertex(GraphicsVector( a, a, a));
        addVertex(GraphicsVector( a, b, a));
        addVertex(GraphicsVector( b, b, a));
        addVertex(GraphicsVector( b, a, a));
        addVertex(GraphicsVector( a, a, b));
        addVertex(GraphicsVector( a, b, b));
        addVertex(GraphicsVector( b, b, b));
        addVertex(GraphicsVector( b, a, b));
 
        // and a list of lines to draw
        static int q[]={0,1,2,3,0,4,5,6,7,4,-1, 1,5,-1,2,6,-1,3,7,-1};
        indexedLineSet(q,sizeof(q)/sizeof(int));
 
        // go to the center and draw the label
        double c = 0.5*(a+b); move_to(GraphicsVector(c,c,c));
        drawText(m_pbox->name());
        set_color("black");
    }
private:
    const SimpleBox* m_pbox;
};
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// implementation of the Demo class

Demo::Demo(gui::DisplayControl & display) : _number(0)
, _box1(new SimpleBox(-400,400, "1"))
, _box2(new SimpleBox(-300,-200, "2"))
, _box3(new SimpleBox(-200,-100, "3"))
, _box4(new SimpleBox(-100, 100, "4"))
{
    using namespace gui;
    
    // create and add their reps to the display, some with buttons in the Display menu
    // to control their visibility
    
    
    DisplayRep* r1 =  new SimpleBoxRep(_box1);
    DisplayRep* r11 = new SimpleBoxRep(_box2);
    DisplayRep* r12 = new SimpleBoxRep(_box3);
    
    DisplayControl::DisplaySubMenu& sm1 = display.subMenu("Boxes..", r1);
    DisplayControl::DisplaySubMenu& sm2 = sm1.subMenu("box 2", r11);
    sm1.add(r12, "box3");
    sm2.add(new SimpleBoxRep(_box4), "box4");
    
    
}

void Demo::process()
{
    // change something that will show up in the display
    ++_number;
}

Demo::~Demo() 
{
    delete _box1;
    delete _box2;
    delete _box3;
    delete _box4;
}

