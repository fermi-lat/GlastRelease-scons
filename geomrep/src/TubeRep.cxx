// $Heading: ConeRep.cxx $
//   Author: T. Burnett
//  Project: Arve graphics

#include "geomrep/TubeRep.h"

#include "geometry/Tube.h"

TubeRep::TubeRep(const Shape& shape)
    : m_tube( dynamic_cast<const Tube& >(shape) )
{};


inline static void createPolyLine(ArvePlottableRep* v, Vector a[], unsigned n)
{
    v->move_to(a[0]); 
    for(unsigned i=1;i<n;i++) 
	v->line_to(a[i]);
}

static inline void createPolyLine(ArvePlottableRep& v, Vector a[], unsigned n)
{
  v.move_to(a[0]); 
  for(unsigned i=1;i<n;i++) 
     v.line_to(a[i]);
}

void 
TubeRep::update()
{
   //  create a representation to draw a Tube by calculating the
   //  vertices of polygons used to display a Tube and storing them
   //  in arrays
   //  the inner and outer cylinders are drawn using ns-sided
   //  polygons at each end, connected by lines from one end to the 
   //  other at the vertices of the polygons
	const int ns = 97, rib=ns/8	;
	Vector vi0[ns];   //  arrays for ns-sided polygons for drawing arcs
	Vector vo0[ns];   //  need 1 set of x, y, z arrays for each end
	Vector vi1[ns];   //  of inner (i) and outer (o) arcs of a Tube
	Vector vo1[ns];
	Vector tv;
	
	int nonzero_inner_radius = m_tube.innerRadius()>0;
//  calculate the coordinates of the inner and outer cylinders
	double hl = 0.5 * m_tube.length();
	double phi = 0.0;
	double dphi = 2.0 * M_PI / ( ns - 1 );
	Vector a = m_tube.axis(),
	       c = m_tube.center(),
	       c0 = c  -  hl * a,
	       c1 = c  +  hl * a;
	Vector c0hat( 1.0, 0.0, 0.0 );
	float ax = a.x();
	float ay = a.y();
	float d = sqrt( ax * ax  +  ay * ay );
	if ( d > 0.0 )
		c0hat = Vector( -ay/d, ax/d, 0.0 );
	Vector c1hat = a.cross( c0hat );
	int i;
	for (  i = 0; i < ns; i++ ) {
		phi = i * dphi;
		tv = cos( phi ) * c0hat  +  sin( phi ) * c1hat;
		vi0[i] = c0  +  m_tube.innerRadius() * tv;
		vi1[i] = c1  +  m_tube.innerRadius() * tv;
		vo0[i] = c0  +  m_tube.outerRadius() * tv;
		vo1[i] = c1  +  m_tube.outerRadius() * tv;
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
	for ( i = 0; i < ns; i += rib ) {
		Vector pl[] = { vo0[i], vo1[i], vi1[i], vi0[i], vo0[i] };
		createPolyLine(this,  pl, 5 );
	}
}

