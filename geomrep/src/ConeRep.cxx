// $Heading: ConeRep.cxx $
//   Author: T. Burnett

#include "geomrep/ConeRep.h"

#include "geometry/Cone.h"


inline static void createPolyLine(gui::DisplayRep* v, Vector a[], unsigned n)
{
    v->move_to(a[0]);
    for(unsigned i=1;i<n;i++)
	v->line_to(a[i]);
}


void ConeRep::update()
{
	const int ns = 25;
	Vector vi0[25];   //  arrays for 24-sided polygons for drawing arcs
	Vector vo0[25];   //  need 1 set of x, y, z arrays for each end
	Vector vi1[25];   //  of inner (i) and outer (o) arcs of a Tube
	Vector vo1[25];
	Vector tv;
	
	int nonzero_inner_radius = m_cone.innerRadius(0)>0;

    //  calculate the coordinates of the inner and outer cylinders
	double hl = 0.5 * m_cone.length();
	double phi = 0.0;
	double dphi = 2.0 * M_PI / ( ns - 1 );

	Vector c = m_cone.center(),
		axis = m_cone.axis(),
	       c0 = c  -  hl * axis,
	       c1 = c  +  hl * axis;
	Vector c0hat( 1.0, 0.0, 0.0 );
	float ax = axis.x();
	float ay = axis.y();
	float d = sqrt( ax * ax  +  ay * ay );
	if ( d > 0.0 )
		c0hat = Vector( -ay/d, ax/d, 0.0 );
	Vector c1hat = axis.cross( c0hat );
	int i;
	for ( i = 0; i < ns; i++ ) {
		phi = i * dphi;
		tv = cos( phi ) * c0hat  +  sin( phi ) * c1hat;
		vi0[i] = c0  +  m_cone.innerRadius(-hl) * tv;
		vi1[i] = c1  +  m_cone.innerRadius(hl) * tv;
		vo0[i] = c0  +  m_cone.outerRadius(-hl) * tv;
		vo1[i] = c1  +  m_cone.outerRadius(hl) * tv;
	}
    //  draw inner cylinder (if inner radius is nonzero)
 	if ( nonzero_inner_radius ) {
		createPolyLine( this, vi0, ns );
		createPolyLine( this, vi1, ns );
	}
    //  draw outer cylinder
	createPolyLine( this, vo0, ns );
	createPolyLine( this, vo1, ns );
   //  draw every third "rib" connecting the polygons
	for ( i = 0; i < ns; i += 3 ) {
		Vector pl[] = { vo0[i], vo1[i], vi1[i], vi0[i], vo0[i] };
		createPolyLine(this,  pl, 5 );
	}

}

