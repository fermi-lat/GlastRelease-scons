// $Header$
//  Author: T. Burnett

#include "geomrep/TrackRep.h"
#include "geometry/Ray.h"

bool TrackRep::useColor=true;
unsigned TrackRep::charged_color_index=1;  //black
unsigned TrackRep::neutral_color_index=0;  //white

float TrackRep::delta = 0.005f;
// maximum error allowed between true trajectory and straight line approx
float TrackRep::minStepSize = 0.5f;


void TrackRep::update()
{
    if( useColor )
        set_col_index( m_track.charged()? charged_color_index : neutral_color_index );
    else
        set_line_style( m_track.charged()? SOLID_LINE : DOTTED_LINE );

    noDetailCheck(); // prevent detail checks
    move_to(m_track.position(0));

    for( Track::const_iterator it= m_track.begin(); it !=m_track.end(); ++it) {
        const Ray* r = *it;
        float stot = r->getArcLength();
        float kappa;
        if( (kappa=r->curvature())!=0 ){
            // helix segment--determine how many steps to take
            float step = sqrt(8*delta/fabs(kappa));  // minimum step for delta
            if( step< minStepSize) step=minStepSize;
            if(stot > 100000) continue;
            int nstep = (int) ( stot / step );
            if(nstep>1) {
                step = stot/nstep;
                float s=0;
                for( int j=0; j<nstep ; j++ )
                    line_to(r->position(s+=step));
            }
        }
        line_to(r->position(stot));   // step to each boundary
    }
    flush();
}


