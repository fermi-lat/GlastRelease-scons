//     $Id$
//  Author: Toby Burnett
//
//

#include "gui/DisplayList.h"

#include "gui/ViewPort.h"
#include "gui/GraphicsVector.h"


#include <cassert>
#include <cfloat>
#include <cmath>

namespace gui {
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  file scope variables
//----------------------
static DisplayList::Node* currentDisp;    // should be either VwPolyLine or VwMarkers
static enum  {NONE, MARKER, POLYLINE, COORDLIST} status ;
static GraphicsVector currentPt;   // current point

//===============================================================================
//                     Limits
//===============================================================================
class DisplayList::Node
// Private abstract class for display primitives
{
friend class DisplayList;
    virtual void flush(){};
    virtual void addPoint(const GraphicsVector&){};
    virtual void draw(ViewPort*)=0;
public:
    virtual ~Node(){};
};

class DisplayList::Limits : public DisplayList::Node {
friend class DisplayList;
    GraphicsVector	limits[2];
    Limits();
    virtual void draw(ViewPort*);
    void update(const GraphicsVector&);
    float checkScale(ViewPort*);
};


DisplayList::Limits::Limits(){
    limits[1] = GraphicsVector(-100,-100,-100);
    limits[0] = GraphicsVector( 100,100,100);
}

static float fmin(float x, float y){return x<y?x:y;}
static float fmax(float x, float y){return x>y?x:y;}

void
DisplayList::Limits::update(const GraphicsVector& pt){
    const GraphicsVector& t = limits[0];     
    const GraphicsVector& u = limits[1];
    limits[0]=GraphicsVector(
		fmin(t.x(), pt.x()),
		fmin(t.y(), pt.y()),
		fmin(t.z(), pt.z()));
    limits[1]=GraphicsVector(
		fmax(u.x(), pt.x()),
		fmax(u.y(), pt.y()),
		fmax(u.z(), pt.z()));
}
void DisplayList::Limits::draw(ViewPort* v){

   v->drawPL(limits,2);
}
float
DisplayList::Limits::checkScale(ViewPort* vw)
{
   float wx1,wx2,wy1,wy2;  // window coords in current view
   if( limits[0].x()==FLT_MAX )
	   return 100;
   vw->transformToWindow(limits[0],wx1,wy1);
   vw->transformToWindow(limits[1],wx2,wy2);
   float dx = wx1-wx2,
       dy = wy1-wy2,
       r = sqrt(dx*dx+dy*dy);
   return r==0? 100 : r;  // exactly zero means a point
}

//===============================================================================

class DisplayList::PolyLine: public DisplayList::Node, public std::vector<GraphicsVector> {
 public:
    PolyLine(const GraphicsVector&);
    virtual ~PolyLine();
    void addPoint(const GraphicsVector&);
    virtual void draw(ViewPort*);
    void flush();
};

//===============================================================================

class DisplayList::Markers : public DisplayList::PolyLine{
  public:
	Markers(const GraphicsVector&);
	virtual void draw(ViewPort*);
};
//===============================================================================

class DisplayList::Text : public DisplayList::Node{
	char * text;
   public:
	Text(const char *);
	virtual void draw(ViewPort*);
	virtual ~Text();
};
//===============================================================================

class DisplayList::Flush : public DisplayList::Node {
   public:
	virtual void draw(ViewPort*);
};
//===============================================================================

class DisplayList::Color : public DisplayList::Node {
	float red,blue,green;
	int colorType;
 public:
	Color(float r,float b,float g,int t){red=r;blue=b;green=g;colorType=t;}
	virtual void draw(ViewPort*);
};
//===============================================================================
class DisplayList::ColIndex : public DisplayList::Node {
        int index;
 public:
        ColIndex(int i){index=i;}
        virtual void draw(ViewPort*);
};
//===============================================================================

class DisplayList::LineWidth :public DisplayList::Node {
	float width;
  public:
	LineWidth(float w){width=w;}
	virtual void draw(ViewPort*);
};
//===============================================================================
class DisplayList::LineStyle : public DisplayList::Node {
  int style;
  public:
	LineStyle(int s){style = s;}
	virtual void draw(ViewPort*);

};
//-----------------------------------------------
//            VRML classes
//-----------------------------------------------

class DisplayList::CoordList : public DisplayList::PolyLine{
  public:
	CoordList(const GraphicsVector&);
	virtual void draw(ViewPort*);
};

//===============================================================================
class DisplayList::Coordinate3 : public DisplayList::Node {
  GraphicsVector* _v;
  int _n;
  public:
      Coordinate3(const GraphicsVector* v , int n);
      ~Coordinate3();
	virtual void draw(ViewPort*);
};
//===============================================================================
class DisplayList::IndexedLineSet : public DisplayList::Node {
 protected:
  int* _ia;
  int _n;
  public:
      IndexedLineSet(const int* ia , int n);
      ~IndexedLineSet();
	virtual void draw(ViewPort*);
};
//===============================================================================
class DisplayList::IndexedFaceSet : public DisplayList::IndexedLineSet {
   public:
      IndexedFaceSet(const int* ia , int n);
	virtual void draw(ViewPort*);
};


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
DisplayList::DisplayList()
//------------------------
   : m_flags(0)
   , m_nodelist(0) , m_nested(0)
{
   currentDisp=0;
   status = NONE;
   currentPt = GraphicsVector(0,0,0);
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
DisplayList::~DisplayList()
//--------------------------
{
    deleteRepresentation();
    if(m_nested) {
	for(DisplayListList::iterator leafptr = m_nested->begin();
	        leafptr != m_nested->end(); ++leafptr)
	    delete *leafptr;
    delete m_nested;
    }
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class DisplayList::Leaf : public DisplayList::Node {
// Special local class to handle draw of nested DisplayList
    // as a primitive
    friend class DisplayList;
    Leaf(DisplayList* dl):m_dl(dl){};
    void draw(ViewPort* v){m_dl->draw(v);}
    DisplayList* m_dl;  
};
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void DisplayList::addDisplayList(DisplayList * v)
//-----------------------------------------------
// put the nested Display list into two lists: the
// m_leaflist for manipulating its visibility, and a special
// DisplayList::Node to draw it in the proper sequence
{
    if( !m_nested )
	m_nested = new DisplayListList;
    m_nested->push_back(v);
    addDisp(new Leaf(v)); //what if no points for size?
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void
DisplayList::append(DisplayList * s)
//----------------------------------
// add the nodes to this one, remove from s
// (ignore any leaves)
{
    if( !s->m_nodelist ) {
	s->update();	// if no representation, have it update
	if( !s->m_nodelist )return; // fail quietly
    }

    DispList::iterator it = s->m_nodelist->begin(),
		end = s->m_nodelist->end();

    // deal with first DisplayList::Node, which sets boundaries
    if( m_nodelist ) { // skip if we have one (should compare)
	delete *it++;  // must delete the Limit object
    }
    else {
	addDisp(0);
	delete (*m_nodelist)[0];
	(*m_nodelist)[0] = *it++;
    }
    // now just copy
    while( it != end )
	addDisp(*it++);

    // now remove them from s, we now own them
    s->m_nodelist->erase(s->m_nodelist->begin(), end);
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void DisplayList::hide(bool outer)
//--------------------------------
{

    if( m_nested	&&  !(*m_nested)[0]->hidden() && ! outer)  {
	for(unsigned i=0; i< m_nested->size(); i++)
	    (*m_nested)[i]->hide();
    }
    else
	m_flags |=  hiddenFlag;
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void DisplayList::show(bool outer)
//--------------------------------
{

    if( hidden() || outer)
	m_flags &= ~hiddenFlag;
    else if( m_nested  )  {
	for(unsigned i=0; i< m_nested->size(); i++)
	    (*m_nested)[i]->show();
    }
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void DisplayList::deleteRepresentation()
//--------------------------------------
{
    if( m_nodelist ){
    //HepAListDeleteAll(*m_nodelist);
	for(unsigned i = 0; i< m_nodelist->size(); ++i)
	    delete (*m_nodelist)[i];
	delete m_nodelist;
    }
    m_nodelist = 0;
    status = NONE;
    if( m_nested)
	for(unsigned i=0; i< m_nested->size(); i++)
	    (*m_nested)[i]->deleteRepresentation();
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void DisplayList::draw(ViewPort* vw)
//------------------------------------------
{
    // check if need to draw
    if( hidden() || (m_nodelist==0) || (m_nodelist->size()==0) ) return;

    DispList::iterator diter = m_nodelist->begin();
    DisplayList::Node* disp = *diter;

    // exit if scale is too small: first DisplayList::Node object does it
    if( ! (m_flags & noCheckDetailFlag)  ){
        float size = ((Limits*)disp)->checkScale(vw);
        if( size < 1 && size < vw->detail() ) {
    	    disp->draw(vw);
	    return;
	}
    }
    if( selected() )  vw->set_enhanced(1);

    // get all the other DisplayList::Node objects to draw themselves

    while( ++diter < m_nodelist->end() )
        (*diter)->draw(vw);

    if( selected() )   vw->set_enhanced(0);

}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void DisplayList::addDisp(DisplayList::Node* d)
//----------------------------------
// internal routine to add a new DisplayList::Node object to the list, and
// do general bookkeeping
// Initialize if m_nodelist is NULL, inserting a Limits
{

    if( !m_nodelist )
    {
	m_nodelist = new DispList; //HepAList<DisplayList::Node>;
	m_nodelist->push_back(new Limits);
    }
    else
	if( currentDisp )
	    currentDisp->flush();

    if( d )
	m_nodelist->push_back(d);
    currentDisp = d;
    status = NONE;

}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void DisplayList::checkScale(const GraphicsVector& pt)
//-----------------------------------------------
{
    DisplayList::Node* first_disp = * (m_nodelist->begin());
    ((Limits*)first_disp)->update(pt);
    currentPt = pt;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void DisplayList::addPoint(const GraphicsVector& pt)
//--------------------------------------------
{
    assert( currentDisp !=0); // FATAL("DisplayList: system error");
    assert( status != NONE) ; // FATAL("DisplayList: no current point");
    currentDisp->addPoint(pt);
    checkScale(pt);
}

//////////////////////////////////////////////////////////////
// Methods that a client can use to create a representation

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void DisplayList::addVertex(const GraphicsVector& pt)
//---------------------------------------------
{
    if( status != COORDLIST){
	addDisp(new CoordList(pt) );
        status = COORDLIST;
    } else addPoint(pt);

}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void DisplayList::coordinate3( const GraphicsVector* v,int n )
//--------------------------------------------------------
{
     addDisp(new Coordinate3(v,n));
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void DisplayList::indexedFaceSet(const int* a, int n)
//---------------------------------------------------
{
     addDisp(new IndexedFaceSet(a,n));
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void DisplayList::indexedLineSet(const int* a, int n)
//---------------------------------------------------
{
     addDisp(new IndexedLineSet(a,n));
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void DisplayList::move_to(const GraphicsVector& pt)
//------------------------------------------
{
        addDisp( new PolyLine(pt) );
        status = POLYLINE;
        checkScale(pt);
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void DisplayList::line_to(const GraphicsVector& pt)
//------------------------------------------
{
    if( status != POLYLINE )
	addDisp( new PolyLine(currentPt) );
    addPoint(pt);
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void DisplayList::markerAt(const GraphicsVector& pt)
//---------------------------------------------
{
    if( status != MARKER){
	addDisp(new Markers(pt) );
        status = MARKER;
        checkScale(pt);
    }
    else
        addPoint(pt);
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void DisplayList::drawText(const std::string& text)
//--------------------------------------------
{
    if( text.size()==0 ) return;
    if( status != POLYLINE )
	addDisp(  new PolyLine(currentPt) );
    addDisp(  new Text(text.c_str()));
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void DisplayList::flush()
//-----------------------
{
    addDisp( new Flush() );
    status = NONE;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void DisplayList::set_rgb(float red, float blue, float green)
//----------------------------------------------------------
{
    addDisp( new Color(red,blue,green,1));
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void DisplayList::set_line_width(float width)
//-----------------------------------------
{
    addDisp( new LineWidth(width) );
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void DisplayList::set_line_style(int style)
{
    addDisp( new LineStyle(style));
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void DisplayList::set_col_index(int i)
//----------------------------------
{
    addDisp( new ColIndex(i));
}

DisplayList::PolyLine::PolyLine(const GraphicsVector& p)
{
   addPoint(p);
}

void
DisplayList::PolyLine::addPoint(const GraphicsVector& pt)
{
   push_back(pt);
}

void
DisplayList::PolyLine::draw(ViewPort* v)
{  // N.B.: this only works for std::vector<GraphicsVector>
    v->drawPL( begin(), size());
}

void
DisplayList::PolyLine::flush()
{
    //pts = (GraphicsVector*)realloc(pts,npts*sizeof(GraphicsVector));
}

DisplayList::PolyLine::~PolyLine()
{
}
//-----------------------------------------------------------------------

DisplayList::Markers::Markers(const GraphicsVector& pt):PolyLine(pt){}

void DisplayList::Markers::draw(ViewPort* v)
{
	v->draw_markers(begin(), size());
}
//-----------------------------------------------------------------------

DisplayList::Text::Text(const char * _text)
{
    text = new char[strlen(_text)+1];
    strcpy( text, _text);
}

void
DisplayList::Text::draw(ViewPort*v)
{
    v->drawText(text);
}

DisplayList::Text::~Text()
{
    delete [] text;

}
//-----------------------------------------------------------------------

void
DisplayList::Flush::draw(ViewPort* v)
{
	v->flush();
}
//-----------------------------------------------------------------------
void
DisplayList::Color::draw(ViewPort* v)
{
       v->set_rgb(red,green,blue);
}
//-----------------------------------------------------------------------
void
DisplayList::ColIndex::draw(ViewPort* v)
{
        v->set_col_index(index);
}
//-----------------------------------------------------------------------
void
DisplayList::LineWidth::draw(ViewPort* v)
{
      v->set_line_width(width);
}
//-----------------------------------------------------------------------
void
DisplayList::LineStyle::draw(ViewPort* v)
{
	v->set_line_style(style);
}
//-----------------------------------------------------------------------
DisplayList::Coordinate3::Coordinate3(const GraphicsVector* v, int n)
: _n(n)
{
   _v = new GraphicsVector[n];
   for(int i=0; i<n; i++)
     _v[i] = v[i]; //note conversion from double to float
}
DisplayList::Coordinate3::~Coordinate3()
{   delete [] _v;
}
void
DisplayList::Coordinate3::draw(ViewPort* v)
{
   v->coordinate3(_v,_n);
}
//-----------------------------------------------------------------------

DisplayList::CoordList::CoordList(const GraphicsVector& pt):PolyLine(pt){}

void
DisplayList::CoordList::draw(ViewPort* v)
{
   v->coordinate3(begin(), size());
}


//-----------------------------------------------------------------------
DisplayList::IndexedLineSet::IndexedLineSet(const int ia[], int n)
:_n(n)
{
   _ia = new int[n];
   for(int i=0; i<n; i++)
     _ia[i]=ia[i];
}
DisplayList::IndexedLineSet::~IndexedLineSet()
{
    delete [] _ia;
}
void
DisplayList::IndexedLineSet::draw(ViewPort* v)
{
    v->indexedLineSet(_ia,_n);
}
//-----------------------------------------------------------------------
DisplayList::IndexedFaceSet::IndexedFaceSet(const int ia[], int n)
 : IndexedLineSet(ia,n){}
void
DisplayList::IndexedFaceSet::draw(ViewPort* v)
{
  v->indexedFaceSet(_ia,_n);
}

}//namespace gui


