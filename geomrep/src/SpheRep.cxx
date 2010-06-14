// $Header$
//  Author: T. Burnett
// Project: Arve graphics
//
// Representation of a Sphe with ArveScene

#include "geomrep/SpheRep.h"
#include "geometry/Sphe.h"

static void createCircle(GraphicsRep& view, const Point& p, double r, const Vector& xhat,
const Vector& yhat)
{
    static const int ns=100;
    static const double dphi= 2.0*M_PI/(ns-1);
    double  phi=0 ;
    view.move_to(p+r*xhat);
    for(int i=0; i<ns; i++)
    {
	phi +=  dphi;
	view.line_to(p+ r*cos(phi)* xhat + r*sin(phi)*yhat);
    }
    view.flush();
}

static const Vector xhat(1,0,0), yhat(0,1,0), zhat(0,0,1);

void SpheRep::update()
{
  double r[]={m_sphe.innerRadius(),m_sphe.outerRadius()};
  for( int i=0; i<2; i++)
  {
    if( r[i]==0 ) continue;
    createCircle(*this, m_sphe.center(), r[i], xhat, yhat );
    createCircle(*this, m_sphe.center(), r[i], yhat, zhat );
    createCircle(*this, m_sphe.center(), r[i], zhat, xhat );
  }
}

