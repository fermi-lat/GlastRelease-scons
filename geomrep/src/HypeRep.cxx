// $Header$
//  Author: Toby Burnett
//
// Display of a Hype with Arve Visualization

#include "geomrep/HypeRep.h"
#include "geometry/Hype.h"


void HypeRep::update()
{
 /*
  //  Create a representation to draw a Hype by calculating the
   //  vertices of polylines used to display a Hype and storing them
   //  in arrays.
   //  To make the drawing less confusing, just draw the inner surface,
   //  with circles represented by 24-sided polygons at each end, connected
   //  by four hyperbolas represented by polylines with 48 segments.
        Vector v1[25];  //  array for 24-sided polygon for +z circle
        Vector v2[25];  //  array for 24-sided polygon for -z circle
        Vector v3[49];  //  array for drawing hyperbolas
//  Calculate vertices of the polygons for the circles in the local coordinate
//  system, then transform them to the global system.
        int i, j, n = 25;
        double phi = 0.0;
        double dphi = 2.0 * M_PI / ( n - 1 );
        double vx, vy;
        Vector vt1, vt2;
        for ( i = 0; i < n; i++ ) {
                phi = i * dphi;
                vx = inner_radius * cos( phi );
                vy = inner_radius * sin( phi );
                vt1 = Vector( vx, vy,  half_length );
                v1[i] = localToGlobal( vt1 );
                vt2 = Vector( vx, vy, -half_length );
                v2[i] = localToGlobal( vt2 );
        }
//  Draw the polylines for the circles.
        createPolyLine( v1, 25 );
        createPolyLine( v2, 25 );
//
//  Now calculate and draw the four hyperbolas.
        Vector vt;
        n = 49;
        double z = -half_length;
        double dz = 2.0 * half_length / n;
        AntiHyperbolic* ah = (AntiHyperbolic *)surface_list->first();
        double ca, sa, rz, lx, ly;
        for ( i = 0; i < 4; i++ ) {
                ca = 0.0;
                sa = 0.0;
                if ( i == 1 )
                        ca = 1.0;
                if ( i == 2 )
                        sa = 1.0;
                if ( i == 3 )
                        ca = -1.0;
                if ( i == 4 )
                        sa = -1.0;
                for ( j = 0; j < n; j++ ) {
                        rz = ah->r_of_z( z );
                        lx = rz * ca;
                        ly = rz * sa;
                        vt = Vector( lx, ly, z );
                        v3[j] = localToGlobal( vt );
                        z += dz;
                }
                createPolyLine( v3, n );
        }

 */
}

