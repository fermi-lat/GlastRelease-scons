// $Header$
//  Author: Hardy, T. Burnett
//
#ifndef WIN32  // stupid to prevent compilation in windows 

#include "Xdraw.h"
#include "gui/DisplayRep.h"
#include <vector>
#include <cmath>

static int markerSize = 4;
static Colormap cmap;
static XColor color_blue, color_red, color_black, colorrgb;
static int useColor=0;  //flag that colors allocated

// table of colors
static std::vector<unsigned> color_table;

int   Xdraw::xwin  (float x)const { return (int)((x+1.)*currentQuad.dx/2) + currentQuad.x; }
int   Xdraw::ywin  (float y)const { return (int)((1.-y)*currentQuad.dy/2) + currentQuad.y; }
float Xdraw::xnorm (int   x)const { return 2 *(x - currentQuad.x) / (float)currentQuad.dx  -  1; }
float Xdraw::ynorm (int   y)const { return 1  -  2 *(y - currentQuad.y) / (float)currentQuad.dy; }

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Xdraw::Xdraw(int width, int height)
  : mydisplay(0)
  , mywindow(0)
{

   // set single/quad rectangles
   flSingle = 1;
   resize (width, height);
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void Xdraw::setDisplay(Display* d, Window w)
{
     mydisplay = d;
     mywindow = w;
     myscreen = DefaultScreen(d);
     if (DisplayPlanes(mydisplay, myscreen) > 1 && !useColor ) {
         useColor = 1;
         cmap = DefaultColormap(mydisplay, myscreen);
         (void)XAllocNamedColor(mydisplay, cmap, "red", &color_red, &colorrgb);
         (void)XAllocNamedColor(mydisplay, cmap, "blue",&color_blue, &colorrgb);
         (void)XAllocNamedColor(mydisplay, cmap, "black",&color_black, &colorrgb);

        // allocate the colors known to GraphicsModel
        // TODO: only define the ones that are used, as requested
        int i=0;
        do{
            XColor color_temp;
            DisplayRep::ColorInfo& cinfo = DisplayRep::pallete[i++];
            if( cinfo.name == 0 ) break;
            (void)XAllocNamedColor(mydisplay, cmap, cinfo.name,
                &color_temp, &colorrgb);
            color_table.push_back( color_temp.pixel );
        }while(1);
    }

}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Xdraw::~Xdraw()
{
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void Xdraw::draw_string(float x, float y, const char *string, int size)
{
    move_to(x,y);
    draw_string(string,size);
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void Xdraw::draw_string(const char* string, int /*size*/)
// draw the string. Note that XDrawImageString makes a light background
{
    XDrawString (mydisplay, mywindow,
		 DefaultGC (mydisplay, myscreen),
		 xwin(lastx), ywin(lasty),
		 string, strlen (string));
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void Xdraw::move_to(float x, float y)
{
    lastx = x;
    lasty = y;
}
static bool inline inside(float x, float y){return fabs(x)<1.01 && fabs(y)<1.01;}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Special clipping routine since X seems to mess up with high magnification
bool Xdraw::clip(float x, float y,
		 int& ix1, int& iy1, int& ix2, int& iy2)
{
  bool b1=inside(x,y), b2=inside(lastx,lasty);
  if( b1 && b2) {
    ix2 = xwin(x);
    iy2 = ywin(y);
    ix1 = xwin(lastx);
    iy1 = ywin(lasty);
    return true; // both endpoints are inside: done
  }
  // either or both is outside
  float dx=lastx-x, dy=lasty-y;
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
  ix1 = xwin(x+t[0]*dx);
  iy1 = ywin(y+t[0]*dy);
  ix2 = xwin(x+t[1]*dx);
  iy2 = ywin(y+t[1]*dy);
  return true;
} 
void Xdraw::line_to(float x, float y)
{
  int x1,y1,x2,y2;
  if( clip(x,y, x1,y1,x2,y2) )
    XDrawLine (mydisplay, mywindow,
	       DefaultGC (mydisplay, myscreen),
	       x1, y1,
	       x2, y2 );
  lastx = x;
  lasty = y;
}
void Xdraw::fill_polygon(const float* xy, int npoints, Shading pattern)
{
  // fill a single polygon
  XPoint * points = new XPoint[npoints];

    XPoint* pt = points;
    for(int i = 0; i<npoints; ++i, ++pt) {
      pt->x = xwin(*xy++);
      pt->y = ywin(*xy++);
    }
    //TODO: figure out how to specify pattern
    XFillPolygon (mydisplay, mywindow, DefaultGC (mydisplay, myscreen),
	       points, npoints, 
	       Convex,CoordModeOrigin);

    delete[] points;

}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void Xdraw::flush()
{

#if 0
   XRectangle rectangle[1];
   rectangle[0].x      = 1;
   rectangle[0].y      = 1;
   rectangle[0].width  = screen.dx-2;
   rectangle[0].height = screen.dy-2;
   XSetClipRectangles(mydisplay,
                      DefaultGC (mydisplay, myscreen),
                      0, 0, rectangle, 1, Unsorted);
#endif
   XFlush(mydisplay);
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void Xdraw::draw_marker(float x, float y)
{
    int thisx = xwin(x);
    int thisy = ywin(y);
 /*
   XSetForeground (mydisplay, DefaultGC(mydisplay,myscreen),
                   color_black.pixel);
*/
    XDrawLine (mydisplay, mywindow,
              DefaultGC (mydisplay, myscreen),
              thisx-markerSize, thisy,
              thisx+markerSize, thisy );
    XDrawLine (mydisplay, mywindow,
              DefaultGC (mydisplay, myscreen),
              thisx, thisy-markerSize,
              thisx, thisy+markerSize );
    XDrawLine (mydisplay, mywindow,
              DefaultGC (mydisplay, myscreen),
              thisx-markerSize, thisy-markerSize,
              thisx+markerSize, thisy+markerSize );
    XDrawLine (mydisplay, mywindow,
              DefaultGC (mydisplay, myscreen),
              thisx-markerSize, thisy+markerSize,
              thisx+markerSize, thisy-markerSize );
/*
  XSetForeground (mydisplay, DefaultGC(mydisplay,myscreen),
                   color_red.pixel);
 */
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void Xdraw::set_quad (int quadnr, const char * title, int /*_selected*/)
{
    if (quadnr==0)      {
       flSingle    = 1;
       currentQuad = screen;
    }   else      {
       flSingle    = 0;
       currentQuad = quadrant [selected = quadnr-1];
    }

    // clear background of view
    XClearArea(mydisplay, mywindow, currentQuad.x, currentQuad.y, currentQuad.dx, currentQuad.dy, False);

    // set clipping for view
    XRectangle rectangle[1];
    rectangle[0].x      = currentQuad.x  + 1;
    rectangle[0].y      = currentQuad.y  + 1;
    rectangle[0].width  = currentQuad.dx - 2;
    rectangle[0].height = currentQuad.dy - 2;
    XSetClipRectangles(mydisplay, DefaultGC (mydisplay, myscreen),
                      0, 0, rectangle, 1, Unsorted);

    // draw window title, if any
    if(title)   XDrawString (mydisplay, mywindow,
                     DefaultGC (mydisplay, myscreen),
                     xwin(-0.9), ywin(-0.9),
                     title, strlen (title));

    set_defaults();
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void Xdraw::set_defaults()
   {
   XSetLineAttributes(mydisplay,
                      DefaultGC (mydisplay, myscreen),
                      0,  // linewidth, Hardy: was former 1
                      LineSolid,
                      CapButt,
                      JoinMiter);
   if (useColor)
      XSetForeground(mydisplay, DefaultGC(mydisplay,myscreen),
                     color_black.pixel);
   }

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void Xdraw::setXorMode(int on)
   {
   if (on)
      {
      /*
      XSetForeground     (mydisplay, DefaultGC(mydisplay,myscreen),
                          color_black.pixel);
                          */
      XSetFunction       (mydisplay, DefaultGC(mydisplay, myscreen),
                          GXxor);
      XSetLineAttributes (mydisplay, DefaultGC(mydisplay, myscreen),
                          0, //1, // linewidth
                          LineSolid,
                          CapButt,
                          JoinMiter);
      }
   else
      {
      /*
      XSetForeground     (mydisplay, DefaultGC(mydisplay,myscreen),
                          color_black.pixel);
                          */
      XSetFunction       (mydisplay, DefaultGC(mydisplay, myscreen),
                          GXcopy);
      XSetLineAttributes (mydisplay, DefaultGC(mydisplay, myscreen),
                          0, //1, // linewidth
                          LineSolid,
                          CapButt,
                          JoinMiter);
      }
   }

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void Xdraw::set_line_style(int style)
   {
   switch (style ){
    case DisplayRep::DOTTED_LINE : style = LineOnOffDash; break;
    case DisplayRep::SOLID_LINE  :
    default:                     style = LineSolid;

   }
  XSetLineAttributes(mydisplay,
                      DefaultGC (mydisplay, myscreen),
                      0, //1, // linewidth
                      style,
                      CapButt,
                      JoinMiter);
   }


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void Xdraw::set_col_index(int index)
{
    XSetForeground(mydisplay, DefaultGC(mydisplay,myscreen), color_table[index]);
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void Xdraw::resize (int xsize, int ysize)
{
    if( xsize==0 )   {
	XWindowAttributes attribs;
	XGetWindowAttributes (mydisplay, mywindow, &attribs);
	xsize = attribs.width; ysize = attribs.height;
    }
    // keep square 
    if( xsize < ysize ) ysize=xsize; else xsize=ysize;
// causes seg. fault??   set_defaults();
   screen.x      = 0;          screen.y      = 0;          screen.dx      = xsize;    screen.dy      = ysize;
   quadrant[0].x = 0;          quadrant[0].y = 0;          quadrant[0].dx = xsize/2;  quadrant[0].dy = ysize/2;
   quadrant[1].x = xsize/2+1;  quadrant[1].y = 0;          quadrant[1].dx = xsize/2;  quadrant[1].dy = ysize/2;
   quadrant[2].x = 0;          quadrant[2].y = ysize/2+1;  quadrant[2].dx = xsize/2;  quadrant[2].dy = ysize/2;
   quadrant[3].x = xsize/2+1;  quadrant[3].y = ysize/2+1;  quadrant[3].dx = xsize/2;  quadrant[3].dy = ysize/2;

   if (flSingle) currentQuad = screen;
   else      currentQuad = quadrant[selected];
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void Xdraw::drawRect(const Xdraw::Rectangle& selection)
{
    move_to( selection.x,                selection.y                );
    line_to( selection.x + selection.dx, selection.y                );
    line_to( selection.x + selection.dx, selection.y + selection.dy );
    line_to( selection.x,                selection.y + selection.dy );
    line_to( selection.x,                selection.y                );
}


#endif //NT_MSVCPP

