// $Header$
//  Author: T. Burnett


#include "Vrml.h"
#include "gui/DisplayRep.h"


namespace gui {

static GraphicsVector currentPoint; // use below to keep track of current point.

static const  char * setup_string = 
"#VRML V1.0 ascii\n"
"DEF BackgroundColor Info {  string  \"0.8 0.8 0.8\"  }\n"
"FontStyle { size 3}\n"
"Material { emissiveColor  0.5  0.5  0.5  ambientColor 0 0 0  }\n";



Vrml::Vrml(std::ostream& c)
:  _out(c)
,  _ncoord(0)
{
    _out << setup_string;
    ViewPort::detail(0);  //should prevent detail check
}
Vrml::~Vrml(){_out.flush();}


//-----------------------------------------------------------------------------
//                     drawing
void Vrml::set_quad(int /*quad*/, int /*selected*/)
{
    
}
void Vrml::beginSeparator(){   _out << "Separator {\n";}
void Vrml::endSeparator()  {   _out << "}\n"; _ncoord=0; }

void Vrml::coordinate3(const GraphicsVector* v,int n)
{
    _out << "Coordinate3 { point [ \n " ;
    for( int i=0; i<n; i++){
        _out << v[i].x()<< ' ' << v[i].y()<< ' '<< v[i].z() ;
        if( i< n-1 ) _out <<",\n"	;
    }
    _out << "] }\n";
    _ncoord=n;
    currentPoint = v[n-1];
}
void Vrml::indexedLineSet(const int ia[], int n)
{
    _out << "IndexedLineSet{ coordIndex[\n";
    for(int i=0; i<n; i++){ 
        _out << ia[i] ;
        if( i< n-1 ) _out << ",";
    }
    _out << "] }\n";
}

void Vrml::indexedFaceSet(const int  *ia,int n)
{
    _out << "IndexedFaceSet{ coordIndex[\n";
    for(int i=0; i<n; i++){ 
        _out << ia[i] ;
        if( i< n-1 ) _out << ",";
    }
    _out << "] }\n";
}
void Vrml::lineSet( int n)
{
    _out << "IndexedLineSet{ coordIndex[\n";
    for(int i=0; i<n; i++) 
        _out << i << ',';
    _out << -1 << "] }\n";
}

void Vrml::drawPL(const GraphicsVector* v,int n)
{  
    coordinate3(v,n);
    lineSet(n);
}

void Vrml::drawText(const char *text, const GraphicsVector& pnt, int /*ht*/)
{
    beginSeparator();
    _out << "Transform { translation "
        << pnt.x() << ' '<< pnt.y() << ' ' << pnt.z() << "}\n";
    _out << "AsciiText { string \"" << text << "\"} \n";
    endSeparator();
}
void Vrml::drawText(const char* text)
{
    drawText(text, currentPoint);
    
}


void Vrml::draw_markers(const GraphicsVector* v,int n)
{ 
    coordinate3(v,n);
    _out << "PointSet { startIndex 0 numPoints " << n << "}\n";
    
}
// draw markers (points) at the specified points

void Vrml::flush(){
    _out.flush();
}
// flush accumulated drawing commands


//-----------------------------------------------------------------------
//  functions for setting drawing attributes

void Vrml::set_rgb( float r, float g, float b )
{
    _out << "Material { emissiveColor "<< r << ' '<< g <<' '<< b << "}\n";
}
void Vrml::set_col_index(int index)
{
    DisplayRep::ColorInfo& cinfo = DisplayRep::pallete[index];
    set_rgb(cinfo.r/255., cinfo.g/255., cinfo.b/255.);
}

void Vrml::set_marker_size(float ){}
void Vrml::set_line_width( float  ){}

void Vrml::set_line_style( int ){}

void Vrml::set_cap_style(int ){}
void Vrml::setJoinStyle(int ){}

void Vrml::setBkColIndex(int ){}
void Vrml::set_enhanced(int){}

//-----------------------------------------------------------------------------
//             stuff associated with projection
/*
void Vrml::setProjection()
{
glMatrixMode(GL_PROJECTION);
float projection[]=
{zoomFactor,0,   0,   0
,0,   zoomFactor*aspectRatio,0,   0
,0,   0,   -zoomFactor, -zoomFactor*perspective
,0,   0,   0,   1};

  glLoadMatrixf(projection);
  
    glMultMatrixf(R);
    
      glScalef(1/scale,1/scale,1/scale);     // overall scale
      glTranslatef(-center->x(),-center->y(),-center->z()); // negative of origin  
      glMatrixMode(GL_MODELVIEW);
      }
*/
} // namespace gui




