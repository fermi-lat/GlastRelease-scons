//     $Id$
//  Author: Toby Burnett
//
// the class DisplayList maintains a list of display commands, with some controls
//
//
#ifndef D0VE_DisplayList_H
#define D0VE_DisplayList_H



#include <string>
#include <vector>
#include "gui/GraphicsVector.h"

namespace gui {

class ViewPort;
class DisplayList;




class DisplayList {

 public:

    DisplayList();
    virtual ~DisplayList();
    // constructor,  destructor

    virtual void draw(ViewPort*);
    // this method may be overridden by specific code

    virtual void update()=0;
    // specific view subclass must implement

    void hide(bool outer=false);
    void show(bool outer=false);
    // turn off/on display m_flags, from inside out if nested
    // if arg is true, do outer only (or all)
   
    void noDetailCheck(){m_flags |= noCheckDetailFlag;} ;
    void select(){  m_flags |=  selectedFlag;}
    void unSelect(){m_flags &= ~selectedFlag;}
    // set/clear hidden, selected, detail checking flag bits

    bool hidden()const{return (m_flags & hiddenFlag)!=0;}
    bool selected()const{return (m_flags & selectedFlag)!=0;}
    // return hidden, selected status


    virtual void deleteRepresentation();
    // delete the representation, if it exists
    // virtual to allow a subclass to modify or bypass

    bool  hasRepresentation()const{return m_nodelist!=0;}
    // report if representation exists

    void append(DisplayList* second);
    // append all the nodes of second to this one.

    void addDisplayList(DisplayList*);
    // add a leaf to this node. Note that if this object is hidden,
    // so are its children


    void changeLevel(int c);
    // use to set level of detail shown; c=0: reset, otherwise change by c, 
    // otherwise stops when gets to the level

    // Following functions save instructions for the default draw

    void moveTo(const GraphicsVector&);
    // move current point to location

    void lineTo(const GraphicsVector&);
    // draw a line to the specified point

    // following are for convenience of those not in a vector world.
    void lineTo(double x, double y, double z){lineTo(gui::GraphicsVector(x,y,z));}
    void moveTo(double x, double y, double z){moveTo(gui::GraphicsVector(x,y,z));}

    void drawText(const std::string&);
    // draw the text at the current point

    void markerAt(const GraphicsVector&);
    // draw a marker at the given point

    void set_rgb(float,float,float);
    // set color according to red, green, blue

    void set_line_width(float);
    // set line width

    void set_line_style(int);
    // set line style (note: need a enum definition of styles)

    void set_col_index(int);
    // set Color index


    void flush();
    // end instructions

    //==== following implement VRML possibility
    void addVertex(const GraphicsVector& v);
    void coordinate3( const GraphicsVector* v,int n );
    // set the Coordinate3 list

    void indexedFaceSet(const int* a, int n);
    void indexedLineSet(const int* a, int n);
    // can draw lines or faces from stored coordinates. See VRML 1.0 spec at
    //   http://vrml.wired.com/vrml.tech/vrml10-3.html

    // nested classes, all inheriting from DisplayNode.  pre declared here
    class Node;
    class Leaf;
    class ColIndex;
    class Color;
    class Coordinate3;
    class CoordList;
    class Flush;
    class IndexedFaceSet;
    class IndexedLineSet;
    class Limits;
    class LineStyle;
    class LineWidth;
    class Markers;
    class PolyLine;
    class Text;

 private:
    typedef std::vector<Node*> DispList ;
    typedef std::vector<DisplayList*> DisplayListList;

    void addDisp(DisplayList::Node*);               // add another displayer
    void addPoint(const GraphicsVector&);
    void checkScale(const GraphicsVector&);

    int     m_flags;          // hidden, checkDetail, selected

    DispList*   m_nodelist;
    // pointer to list of pointer to  primitive DisplayNode objects

    DisplayListList * m_nested;
    // pointer to list of internal DisplayList

    enum Flags {hiddenFlag = 1,  selectedFlag = 2, noCheckDetailFlag=4};
    // used to interpret bits in the m_flags word
};

} // namespace gui
#endif
