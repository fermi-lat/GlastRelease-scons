//  $Header$
//   Author: T. Burnett

#include "geomrep/ConsRep.h"

#include "geometry/Cons.h"


inline static void createPolyLine(gui::DisplayRep* v, Vector a[], unsigned n)
{
    v->move_to(a[0]);
    for(unsigned i=1;i<n;i++)
	v->line_to(a[i]);
}


void ConsRep::update()
{

}

