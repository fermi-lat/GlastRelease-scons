//     $Id$
//  Author: Toby Burnett
//
//
// THis is a ViewPort subclass that writes 3-d information in VRML format to a file.
// The file can then be viewed with various WWW browsers.
//
// sample usage:
//
//   Vrml scene("myscene.wrl"); // will write to the file
//   DisplayList* detector;        //from somewhere
//   detector->draw(&scene);
//
#ifndef VRML_H
#define VRML_H

#include "gui/ViewPort.h"


#include <string>
#include <iostream>

namespace gui {
    class Vrml  : public ViewPort
    {
    public:
        Vrml(std::ostream& out);
        virtual ~Vrml();
        
        //  functions for setting drawing attributes
        
        void set_rgb( float red, float green, float blue );
        void set_line_width( float lwid );
        void set_line_style( int s );
        
        void set_col_index(int index);
        void set_cap_style(int style);
        void setJoinStyle(int style);
        void setBkColIndex(int index);
        void set_marker_size(float size);
        void set_enhanced(int);
        // set the selected attribute on or off
        
        // draw commands
        
        void set_quad(int quadrant, int highlighted);
        
        void drawPL(const GraphicsVector*,int);
        // draw a polyline 
        
        void drawText(const char *text, const GraphicsVector& pnt, int ht=0);
        void drawText(const char* text);
        // draw (horizontal) text at the given, or current point
        
        void draw_markers(const GraphicsVector*,int);
        // draw markers (points) at the specified points
        
        void beginSeparator();
        void endSeparator();
        // useful for VRML color, texture, etc.
        
        void coordinate3(const GraphicsVector *points, int n);
        // set the points to display
        
        void indexedLineSet(const int* a, int n);
        void indexedFaceSet(const int* a, int n);
        // create polylines or set of faces from coordinate3 list
        
        void flush();
        // flush accumulated drawing commands
        
        
    private:
        
        void lineSet(int n);
        // create a single polyline from all the coordinates
        
        std::ostream& _out;
        int _ncoord; 
        // number of valid coords
        
    };
    
} //namespace gui

#endif 

