//  $Header$
// PostScript drawer imlementation 
// created by Hardy, 28 Mar 94

#include "PSdraw.h"
#include "gui/DisplayRep.h"


#include <iostream>
#include <iomanip>

namespace gui {

static int markerSize = 10;


inline int PSdraw::xwin(float x)
{
    return (int)((x+1)*currentQuad.dx/2) + currentQuad.x;
}

inline int PSdraw::ywin(float y)
{
    return (int)((y+1)*currentQuad.dy/2) + currentQuad.y;
}


struct PSdraw::PageFormat
{
    const char *format;
    int   ptXsize;
    int   ptYsize;
    
    static const PageFormat *liste;
    static void lookUp (const char *_format, int *_ptXsize, int *_ptYsize);
};


void PSdraw::PageFormat::lookUp (const char *_format, int *_ptXsize, int *_ptYsize)
{
    
    for (int i=0; liste[i].format; i++)
    {
        
        if (strcmp(liste[i].format, _format)==0)
        {
            *_ptXsize = liste[i].ptXsize;
            *_ptYsize = liste[i].ptYsize;
            return;
        }
    }
    *_ptXsize = liste[0].ptXsize;
    *_ptYsize = liste[0].ptYsize;
    return;
}


static const PSdraw::PageFormat formatListe[] = 
{
    { "DinA4",       2100, 2970 },
    { "DinA4-quer",  2970, 2100 },
    { 0,             0,    0    }
};
const PSdraw::PageFormat *PSdraw::PageFormat::liste = formatListe;


static const char *psInitString = 
"%!PS-Adobe-2.0\n"
"%%Title:        gui-Display (Portrait A 4)\n"
"%%Pages:        atend\n"
"%%Creator:      Hardy's PS-Printer 08.15\n"
"%%CreationDate: 94/03/28   10.13\n"
"%%EndComments\n"
"/s {stroke} def /l {lineto} def /m {moveto} def /t {translate} def\n"
"/sw {stringwidth} def /r {rlineto} def /rl {roll} def\n"
"/gr {grestore} def\n"
"/c {setrgbcolor} def /lw {setlinewidth} def /sd {setdash} def\n"
"/cl {closepath} def /sf {scalefont setfont} def\n"
"/rgb {setrgbcolor} def\n"
"\n"
"/Helvetica findfont\n"
"30 scalefont\n"               // THB: was 10
"setfont\n"
"0.28346456 0.28346456 scale\n"           // change to 1/10 mm
"1 setlinewidth\n"
"gsave\n";


static const char *psFlushString = 
"";

static const char *psCloseString = 
"\ngrestore\n"
"showpage gr\n"
"%%Trailer\n"
"%%Pages: 1\n"
"gr gr\n"
"%%EOF\n";




PSdraw::PSdraw (const char *filename, const char *page_format)
: psfile( filename )
{
    // send initstring
    psfile << psInitString;
    
    
    // get format for mapping
    PageFormat::lookUp (page_format, &screenWidth, &screenHeight);
    
    
    // preserve aspect ratio
    flSingle = 1;
    resize(screenWidth, screenHeight);
    
    // no stroke necessary, yest
    isStrokePending = 0;
}


PSdraw::~PSdraw()
{
    flush ();
    psfile << psCloseString;
    psfile.flush();
}


void PSdraw::draw_string (float x, float y, const char *string, int size)
{
    move_to (x,y);
    draw_string (string, size);
}


void PSdraw::draw_string (const char* string, int /*size*/)
{
    psfile << "(" << string << ") show\n";
}


void PSdraw::move_to (float x, float y)
{
    lastx = xwin(x); 
    lasty = ywin(y);
    
    if (isStrokePending)
    {
        psfile << " s\n";
        isStrokePending = 0;
    }
    
    psfile << " " << lastx << " " << lasty << " m";
}


void PSdraw::line_to (float x, float y)
{
    int thisx = xwin(x);
    int thisy = ywin(y);
    
    psfile << " " << thisx << " " << thisy << " l";
    isStrokePending = 1;
    
    lastx = thisx;
    lasty = thisy;
}


void PSdraw::flush ()
{
    if (isStrokePending)
    {
        psfile << " s\n";
        isStrokePending = 0;
    }
    
    psfile << psFlushString;
    psfile.flush ();
}


void PSdraw::draw_marker (float x, float y)
{
    if (isStrokePending)
    {
        psfile << " s\n";
        isStrokePending = 0;
    }
    
    int thisx = xwin(x);
    int thisy = ywin(y);
    
    psfile << " " << thisx-markerSize << " " << thisy+markerSize << " m";
    psfile << " " << thisx+markerSize << " " << thisy-markerSize << " l";
    psfile << " " << thisx-markerSize << " " << thisy-markerSize << " m";
    psfile << " " << thisx+markerSize << " " << thisy+markerSize << " l";
    psfile << " " << lastx << " " << lasty << " m";
    
}


void PSdraw::set_quad (int quadnr, const char * title, int /*_selected*/)
{
    if (quadnr==0)
    {
        flSingle    = 1;
        currentQuad = screen;
    }
    else
    {
        flSingle    = 0;
        currentQuad = quadrant [selected = quadnr-1];
    }
    
    // clear background of view ?
    
    // show clipping rectangle
    psfile << "\ngrestore "
        << currentQuad.x   << ' ' << currentQuad.y  << " m "
        << currentQuad.dx  << ' ' << 0              << " r "
        << 0               << ' ' << currentQuad.dy << " r "
        << -currentQuad.dx << ' ' << 0              << " r "
        << " closepath stroke\n";
    
    // set clipping for view
    psfile << "gsave "
        << currentQuad.x   << ' ' << currentQuad.y  << " m "
        << currentQuad.dx  << ' ' << 0              << " r "
        << 0               << ' ' << currentQuad.dy << " r "
        << -currentQuad.dx << ' ' << 0              << " r "
        << " closepath clip\n";
    // enter title
    psfile << xwin(-0.9f) << ' ' << ywin(-0.9f) << " m\n"; // rmoveto\n";
    psfile << '('<<title << ") show\n"; // NOTE: should scan for parens
    
    
}


void PSdraw::resize (int xsize, int ysize)
{
    int xoffset = 0;
    int yoffset = 0;
    
    if (ysize > xsize)
    {
        yoffset = (ysize-xsize) / 2;
        ysize   = xsize;
    }
    else
    {
        xoffset = (xsize-ysize) / 2;
        xsize   = ysize;
    }
    
    screen.x      = 0;          screen.y      = 0;          screen.dx      = xsize;    screen.dy      = ysize;
    quadrant[0].x = 0;          quadrant[0].y = 0;          quadrant[0].dx = xsize/2;  quadrant[0].dy = ysize/2;
    quadrant[1].x = xsize/2+1;  quadrant[1].y = 0;          quadrant[1].dx = xsize/2;  quadrant[1].dy = ysize/2;
    quadrant[2].x = 0;          quadrant[2].y = ysize/2+1;  quadrant[2].dx = xsize/2;  quadrant[2].dy = ysize/2;
    quadrant[3].x = xsize/2+1;  quadrant[3].y = ysize/2+1;  quadrant[3].dx = xsize/2;  quadrant[3].dy = ysize/2;
    
    screen.x      += xoffset;    screen.y      += yoffset;
    quadrant[0].x += xoffset;    quadrant[0].y += yoffset;
    quadrant[1].x += xoffset;    quadrant[1].y += yoffset;
    quadrant[2].x += xoffset;    quadrant[2].y += yoffset;
    quadrant[3].x += xoffset;    quadrant[3].y += yoffset;
    
    
    if (flSingle)
        currentQuad = screen;
    else
        currentQuad = quadrant[selected];
}


void PSdraw::set_defaults ()
{
    //XSetLineAttributes(mydisplay,
    //	      DefaultGC (mydisplay, myscreen),
    //	      1,
    //	      LineSolid,
    //	      CapButt,
    //	      JoinMiter);
    
    //if (useColor) 
    //   XSetForeground(mydisplay, DefaultGC(mydisplay,myscreen), color_black.pixel);
}


void PSdraw::set_line_style (int /*style*/)
{
    //XSetLineAttributes(mydisplay,
    //	      DefaultGC (mydisplay, myscreen),
    //	      1,  // linewidth
    //	      style,
    //	      CapButt,
    //	      JoinMiter);
}
//----------------------------------------------------------------
//                         color
//----------------------------------------------------------------
void PSdraw::set_rgb(float r,float g,float b)
{
  psfile << std::setprecision(3) 
	 << ' ' << r/255. << ' ' << g/255. << ' ' << b/255. << " rgb\n";
}
void PSdraw::set_col_index(int i)
{
    DisplayRep::ColorInfo& cinfo = DisplayRep::pallete[i];
    set_rgb(cinfo.r, cinfo.g, cinfo.b );
    
}



}//namespace gui
