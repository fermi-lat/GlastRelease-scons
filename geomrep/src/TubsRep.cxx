// $Header$
//  Author: Toby Burnett
//
// Display of a Tubs

#include "geomrep/TubsRep.h"
#include "geometry/Tubs.h"


static void createPolyLine(gui::DisplayRep* v, Vector a[], unsigned n)
{
  v->move_to(a[0]);
  for(unsigned i=1;i<n;i++)
     v->line_to(a[i]);
}

void TubsRep::update()
{
    // We need two vectors that will point the direction that the two wedge pieces
    // are going away from the axis.
    // Note that the axis of the wedge is always the axis of the cylinder
    Vector axis = m_tubs.axis(),
           p1 = - m_tubs.wedge().n1().cross(axis),
           p2 =   m_tubs.wedge().n2().cross(axis);

    double openAngle = acos(p1*p2);
    if( (p1.cross(p2))*axis<0 )
	openAngle += M_PI;

    // Let's set up some arrays to hold the points in...
    const int ns = 25; // Number of sides to the circly part
    Vector inside1[ns], inside2[ns];
    Vector outside1[ns + 2], outside2[ns +2];
    // need two extra here for the first and last lines, which will connect the inside and
    // outside parts

    // Let's begin going around the circle.....
    double dphi = 2.0*M_PI/ns;
    double in = m_tubs.innerRadius();
    in = ((in>0)? in:0);
    double out = m_tubs.outerRadius();

    Vector dir = p1,
           disp = axis*.5*m_tubs.length(),
           center = m_tubs.center();


    int c = -1, i; // c is the counter

    for( i=0;  i< openAngle/dphi; i++ ) {
        c++;
        inside1[c] = center + in*dir + disp;
        inside2[c] = center + in*dir - disp;
        outside1[c+1] = center + out*dir + disp;
        outside2[c+1] = center + out*dir - disp;
	dir.rotate(dphi,axis);
    }

    // Draw the last segment so it is complete
    c++;
    dir = p2;
    inside1[c] = center + in*dir + disp;
    inside2[c] = center + in*dir - disp;
    outside1[c+1] = center + out*dir + disp;
    outside2[c+1] = center + out*dir - disp;


    // Now that we have the edges, let's put the two sides in
    outside1[0] = inside1[0];
    outside1[c+2] = inside1[c];
    outside2[0] = inside2[0];
    outside2[c+2] = inside2[c];

    // Now draw the polygon
    createPolyLine(this, outside1, c+3);
    createPolyLine(this, inside1, c+1);
    createPolyLine(this, outside2, c+3);
    createPolyLine(this, inside2, c+1);

    // Now we need some connecting points to make it have depth....
    Vector line_seg[2];
    for ( i=0;i<c;i+=3)
    {
       line_seg[0] = outside1[i+1];
       line_seg[1] = outside2[i+1];
       createPolyLine(this, line_seg, 2);
       line_seg[0] = inside1[i];
       line_seg[1] = inside2[i];
       createPolyLine(this, line_seg, 2);
    }

    // Now we need to connect the end
    line_seg[0] = outside1[c+1];
    line_seg[1] = outside2[c+1];
    createPolyLine(this, line_seg, 2);

    line_seg[0] = inside1[c];
    line_seg[1] = inside2[c];
    createPolyLine(this, line_seg, 2);


}

