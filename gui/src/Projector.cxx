//     $Id$
//  Author: Toby Burnett
//

#ifdef _MSC_VER
#pragma warning(disable: 4244)
#endif

#include "gui/Projector.h"
#include "gui/GraphicsVector.h"

namespace gui {

Projector::Projector(Draw2D*d, const GraphicsVector* ref, float sc,
		     const float * eulers, const char* t)
		     : ViewPort(ref, sc, eulers, t)
		     , m_canvas(d)
{
}
// constructor

Projector::Projector(Draw2D *d)
: m_canvas(d)
{ 
}

ViewPort& 
Projector::copy(const GraphicsVector* ref, float sc,
		const float * eulers, const char* t)
{
    Projector& cpy = *new Projector(m_canvas, ref, sc, eulers, t);
    return cpy;
}

void 
Projector::set_quad(int quadrant, int highlighted)
{      m_canvas->set_quad(quadrant, title(), highlighted);
}

void 
Projector::drawText(const char *text, const GraphicsVector& point, int ht)
{
    float   ox, oy;
    transformToWindow(point, ox, oy);
    m_canvas->draw_string( ox, oy, text, ht);
}


void Projector::draw_markers(const GraphicsVector* pts,int ln)
{
    for( int i=0; i<ln; i++){
	float winx, winy;
	transformToWindow(pts[i],winx, winy);
	m_canvas->draw_marker(winx,winy);
    }
    
}
void Projector::drawPL(const GraphicsVector* pts, int ln)
{
    float winx, winy;
    
    m_canvas->new_path();
    
    transformToWindow(*pts, winx, winy);
    m_canvas->move_to(winx,winy);
    
    for( int i=1; i<ln; i++){
	transformToWindow(pts[i], winx, winy);
	m_canvas->line_to(winx, winy);
    }
    m_canvas->stroke();
    
}
void Projector::drawText(const char* text)
{
    m_canvas->draw_string(text);
}

//-------------------------------------------------------------------
class Point2D {
public:
    float x,y;
};   

static Point2D* pt2=0;

void Projector::coordinate3(const GraphicsVector *pt3, int n)
{
    if( pt2 )
	delete [] pt2;
    pt2 = new Point2D[n];
    for(int i=0; i<n; i++)
	transformToWindow(pt3[i], pt2[i].x, pt2[i].y);
    
}

void Projector::indexedFaceSet(const int* a, int n)
{
    int newseg=1;
    int st=-1;
    while(n--)
    {	int i = *a++;
    if( newseg){
	if( st >=0 )
	    m_canvas->line_to(pt2[st].x,pt2[st].y);
	m_canvas->move_to(pt2[i].x, pt2[i].y);
	st=i;
	newseg=0;
    }
    else if( i>= 0 )
	m_canvas->line_to(pt2[i].x, pt2[i].y);
    else
	newseg=1;
    }
    if( st >=0 )
	m_canvas->line_to(pt2[st].x,pt2[st].y);
    
}
void Projector::indexedLineSet(const int* a, int n)
{
    int newseg=1;
    while(n--)
    {	int i = *a++;
    if( newseg){
	m_canvas->move_to(pt2[i].x, pt2[i].y);
	newseg=0;
    }
    else if( i>= 0 )
	m_canvas->line_to(pt2[i].x, pt2[i].y);
    else
	newseg=1;
    }
}
void Projector::flush()
{
    m_canvas->flush();
    if(pt2)
	delete [] pt2;
    pt2=0;
}

} // namespace gui
