//  $Header$
//   Author: Toby Burnett
//
#ifdef WIN32  // stupid to prevent compilation on unix platforms 
#include "WinDraw.h"
#include "gui/DisplayRep.h"
#include <cassert>
#include <vector>
#include <cmath>

static float scale;

// specific brushes should be done more generally
static HBRUSH hbrBlack=0, hbrRed=0, hbrGreen=0, 
hbrBlue=0, hbrDashed, hbrDotted, hbrLtGray;
static HGDIOBJ hpenOld=0,hpen;

// table of brushes created from pallete table
static std::vector<HBRUSH > hbrTable;

// rectangles used in display
static RECT screen;
static RECT currentQuad;
static RECT quadrant[4];
static RECT selectedQuad;
static bool flSingle=true;

inline int xscale(float x){ return (int)((x+1.0)*scale)+currentQuad.left;}
inline int yscale(float y){ return (int)((1.0-y)*scale)+currentQuad.top;}

static bool inline inside(float x, float y){return fabs(x)<1.01 && fabs(y)<1.01;}
// Special clipping routine to fix high magnification
bool WinDraw::clip(float x, float y,
		 int& ix1, int& iy1, int& ix2, int& iy2)
{
  bool b1=inside(x,y), b2=inside(m_lastx,m_lasty);
  if( b1 && b2) {
    ix2 = xscale(x);
    iy2 = yscale(y);
    ix1 = xscale(m_lastx);
    iy1 = yscale(m_lasty);
    return true; // both endpoints are inside: done
  }
  // either or both is outside
  float dx=m_lastx-x, dy=m_lasty-y;
  // find parametric values for intersections with scene boundaries at +-1
  float s[]={-1,-1,-1,-1};
  int k=0; // count of intersections
  if( dx!=0){s[k++]=(1-x)/dx; s[k++] = (-1-x)/dx;}
  if( dy!=0){s[k++]=(1-y)/dy; s[k++] = (-1-y)/dy;}

  float t[4];
  int j=0;
  if( b1) t[j++]=0;
  if( b2) t[j++]=1;
  
  // if endpoint in region, already have it. now find 1 or 2 others
  for( int i=0; i<k; ++i){
    if(s[i]>0 && s[i]<1 && inside(x+s[i]*dx,y+s[i]*dy)) t[j++]=s[i];
  }
  if( j<2 )return false;
  ix1 = xscale(x+t[0]*dx);
  iy1 = yscale(y+t[0]*dy);
  ix2 = xscale(x+t[1]*dx);
  iy2 = yscale(y+t[1]*dy);
  return true;
} 

void WinDraw::move_to(float x, float y)
{
  m_lastx=x; m_lasty=y;
}

void WinDraw::line_to(float x, float y)
{
  int x1,y1,x2,y2;
  if( clip(x,y, x1,y1,x2,y2) ) {
      ::MoveToEx(m_hdc, x1,y1, NULL);
      ::LineTo(m_hdc, x2,y2);
  }   
  m_lastx = x;
  m_lasty = y;

}

float WinDraw::xnorm(int ix)
{
    return (ix-selectedQuad.left)/scale-1.0;
}
float WinDraw::ynorm(int iy)
{
    return 1.0 - (iy-selectedQuad.top)/scale;
}

WinDraw::WinDraw(HWND win):m_win(win){};
WinDraw::~WinDraw() {}


static float fmin(float x, float y){return x<y?x:y;}
static float fmax(float x, float y){return x>y?x:y;}

void
WinDraw::setHDC(HDC hDC)
{   m_hdc = hDC;
    ::GetClientRect(m_win, &screen);

    int xsize = screen.right, 
	ysize = screen.bottom;
    quadrant[0].left = 0;          quadrant[0].top = 0;          quadrant[0].right = xsize/2;  quadrant[0].bottom = ysize/2;
    quadrant[1].left = xsize/2+1;  quadrant[1].top = 0;          quadrant[1].right = xsize;    quadrant[1].bottom = ysize/2;
    quadrant[2].left = 0;          quadrant[2].top = ysize/2+1;  quadrant[2].right = xsize/2;  quadrant[2].bottom = ysize;
    quadrant[3].left = xsize/2+1;  quadrant[3].top = ysize/2+1;  quadrant[3].right = xsize;    quadrant[3].bottom = ysize;
    
    currentQuad = screen;
    scale = 0.5*fmin(xsize, ysize);
    flSingle = true;
    
    if( hbrBlack ) return;
    
    // define all pens using the static pallete defined in DisplayRep
    // TODO: only create pens as needed
    int i=0;
    do{ 
	DisplayRep::ColorInfo& cinfo = DisplayRep::pallete[i++];
	if( cinfo.name == 0 ) break;
	hbrTable.push_back(
	    (HBRUSH)CreatePen( PS_SOLID, 1,RGB(cinfo.r,cinfo.g,cinfo.b) ) 
	    );
    }while(1);
    
    // (old style pens)
    hbrBlack = (HBRUSH)GetStockObject(BLACK_PEN);
    hbrLtGray =(HBRUSH)GetStockObject (LTGRAY_BRUSH);
    hbrDashed = (HBRUSH)CreatePen(PS_DASH, 1, RGB(0,0,0));
    hbrDotted = (HBRUSH)CreatePen(PS_DOT, 1, RGB(0,0,0));
    hbrRed =   (HBRUSH)CreatePen(PS_SOLID, 2, RGB(255,0,0));
    hbrGreen = (HBRUSH)CreatePen(PS_SOLID, 2, RGB(0,255,0));
    hbrBlue =  (HBRUSH)CreatePen(PS_SOLID, 2, RGB(0,0,255));
}

void WinDraw::drawTitle(const char * title, int selected)
{
    //draw_string(0.f,0.95f,title);
    draw_string(-0.90f,-0.90f,title);
}
void WinDraw::set_quad(int quadnr, const char* title, int selected)
{
    // select the specified quadrant to draw in
    if (quadnr==0)      {
	
	// 0 -- full screen with one view
	flSingle    = true;
	currentQuad = screen;
	selectedQuad = screen;
    }else {
	// 1-4-- drawing view quadnr in corresponding screen quadrant
	flSingle    = false;
	currentQuad = quadrant [ quadnr-1];
	if( selected) selectedQuad = currentQuad;
    }
    
    
    // set to clip to the dimensions of the region
    ::BeginPath(m_hdc);
    ::Rectangle(m_hdc,currentQuad.left, currentQuad.top, currentQuad.right, currentQuad.bottom);
    ::EndPath(m_hdc);
    ::SelectClipPath(m_hdc, RGN_COPY); 
    
    
    scale = 0.5*fmin(currentQuad.right-currentQuad.left, currentQuad.bottom-currentQuad.top);
    
    // if title specfied, draw  the title
    if( title )
        drawTitle(title, selected);
    if(selected) {
	drawRectangle();  
    }
}

void WinDraw::drawRectangle()
{
    static int b=1, c=2; // pixels in from borders
    ::SelectObject(m_hdc, hbrBlack);
    ::MoveToEx(m_hdc,currentQuad.left+b,    currentQuad.top+b,	 NULL);
    ::LineTo(m_hdc,  currentQuad.right-c,    currentQuad.top+b);
    ::LineTo(m_hdc,  currentQuad.right-c,    currentQuad.bottom-c);
    ::LineTo(m_hdc,  currentQuad.left+b,	    currentQuad.bottom-c);
    ::LineTo(m_hdc,  currentQuad.left+b,	    currentQuad.top+b);
    
}
void WinDraw::fill_polygon(const float* xy, int npoints, Shading pattern)
{
    tagPOINT* points = new tagPOINT[npoints];
    tagPOINT* pt = points;
    for(int i = 0; i<npoints; ++i, ++pt) {
	pt->x = xscale(*xy++);
	pt->y = yscale(*xy++);
    }
    if( pattern==NONE){
	::MoveToEx(m_hdc, points->x, points->y,NULL);
	::PolylineTo(m_hdc, points, npoints);
	::LineTo(m_hdc, points->x, points->y);
    }else {
	switch (pattern) {
	case BRIGHT: ::SelectObject(m_hdc, (HBRUSH)GetStockObject (WHITE_BRUSH)); break;
	case MEDIUM: ::SelectObject(m_hdc, hbrLtGray); break;
	case DARK:   ::SelectObject(m_hdc, (HBRUSH)GetStockObject (BLACK_BRUSH));
	}
	::Polygon(m_hdc, points, npoints);
    }
    delete points;

}

void WinDraw::clearArea()
{
}


void WinDraw::draw_string(const char* s, int)
{ 
    const int vertOffset=-5; // arbitray to center vertically
    // ::TextOut(m_hdc,lastx,lasty+vertOffset,s,strlen(s));
    HFONT hfnt= (HFONT)GetStockObject(ANSI_VAR_FONT),
	hOldFont = (HFONT)SelectObject(m_hdc, hfnt);
    COLORREF hOldbk = (COLORREF)::SetBkColor(m_hdc, RGB(192,192,192) );

    //uncomment for center ::SetTextAlign(m_hdc,TA_CENTER ); 
    ::TextOut(m_hdc,xscale(m_lastx), yscale(m_lasty)+vertOffset,s,strlen(s));
    SelectObject(m_hdc, hOldFont); 
    
}
void WinDraw::draw_string(float x,float y,const char *s, int)
{
    move_to(x,y);
    draw_string(s);
}

void WinDraw::clearArea(int x, int y, int width, int height)
{
}


static float ms=0.015f;
void WinDraw::draw_marker(float x, float y)
{
    move_to(x-ms,y);    line_to(x+ms,y);
    move_to(x,   y-ms); line_to(x,   y+ms);
    move_to(x-ms,y-ms); line_to(x+ms,y+ms);
    move_to(x+ms,y-ms); line_to(x-ms,y+ms);
}


void WinDraw::set_defaults()
{   
    ::SelectObject(m_hdc, hbrBlack);
}

void WinDraw::set_line_width(float w)
{
}

// set the line width
void WinDraw::set_line_style( int s)
{
    switch( s)
    {case 0/*SOLID_LINE*/:  
    ::SelectObject(m_hdc, hbrBlack);
    break;
    case 1/*DOTTED_LINE*/: 
	::SelectObject(m_hdc, hbrDotted);
	break;
    case 3/*DASHED_LINE*/: 
	::SelectObject(m_hdc, hbrDashed);
	break;
    case 2/*CENTER_LINE*/:  break;
    }
}

void WinDraw::set_col_index(int i)
{
    assert( i<hbrTable.size() );
    ::SelectObject(m_hdc, hbrTable[i] );
}

void WinDraw::set_rgb(float r, float g, float b)
{
    if( r==0 && g==0 && b==0 )	       hpen = hbrBlack;
    else if (r==1 && g==0 && b==0)      hpen = hbrRed;
    else if (r==0 && g==1 && b==0)      hpen = hbrGreen;
    else if (r==0 && g==0 && b==1)      hpen = hbrBlue;
    else if (r==1 && g==1 && b==1)      hpen = hbrBlack;  // why here?
    else
	hpen =  CreatePen(PS_SOLID, 1, 
	RGB(static_cast<int>(r*256.),
	static_cast<int>(g*256.),
	static_cast<int>(b*256.)));
    
    ::SelectObject(m_hdc, hpen);
    
}

void WinDraw::set_marker_size(float ){
}

void WinDraw::invalidate()
{
    ::InvalidateRect(m_win, NULL, TRUE);
}


#endif // NT_MSVCPP
