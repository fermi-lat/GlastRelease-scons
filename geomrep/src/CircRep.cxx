// $Header$
//

#include "geomrep/CircRep.h"

/*static */void CircRep::update()
{
gui::DisplayRep& view=*this;
    static const int ns=100;
    static const double dphi= 2.0*M_PI/(ns-1);
    double  phi=0 ;
    view.move_to(m_p+m_r*m_xhat);
    for(int i=0; i<ns; i++)
    {
	phi +=  dphi;
	view.line_to(m_p+ m_r*cos(phi)* m_xhat + m_r*sin(phi)*m_yhat);
    }
    view.flush();
}

