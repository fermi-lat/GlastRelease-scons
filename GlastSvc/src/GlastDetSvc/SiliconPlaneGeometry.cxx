//$Header$
#include "SiliconPlaneGeometry.h"

#include <cmath>
#include <algorithm>
/*
In xmlGeoDbs/xml/flight/flightTKRCountPrim.xml
   
      nFeChips = 6
      stripPerWafer = 384
      nWaferAcross = 4

In xmlGeoDbs/xml/flight/flightTKRDimPrim.xml
  
      ssdGap = 0.025
      SiWaferSide = 89.5
      SiWaferActive = 87.552
      ladderGap = 0.20
Will assume that 
  guard ring = 0.5*(SiWaferSide-SiWaferActive) 
*/

// static variable implementations
unsigned int SiliconPlaneGeometry::s_n_si_dies=4; // # of Si dies across single SSD plane
unsigned int SiliconPlaneGeometry::s_stripPerWafer =384;
double       SiliconPlaneGeometry::s_dice_width     =89.500; // width of a single silicon die 

// width of the dead area around the perimeter of a die (from edge to active)
double       SiliconPlaneGeometry::s_guard_ring     = 0.974;

// width of the gap between SSD dies (glue & wirebonds, etc.) (cm)
double       SiliconPlaneGeometry::s_ssd_gap        = 0.025;            
double       SiliconPlaneGeometry::s_ladder_gap     = 0.200;


// insideActiveArea - does a point fall on (in local coords) an active 
// part of the detector? (-) if no
double SiliconPlaneGeometry::insideActiveArea (double x, double y) 
{
    double insideX = insideActiveArea1D(x, ladder_gap());
    double insideY = insideActiveArea1D(y, ssd_gap());
    
    if (insideX>0 || insideY>0)  {
        return std::min(insideX, insideY);
    } else {
        // special case, want distance to the nearest corner?
        // probably doesn't matter much
        return -sqrt(insideX*insideX + insideY*insideY);
    }		
}


/// number of silicon dies
double   SiliconPlaneGeometry::die_width () {     return s_dice_width;}
double   SiliconPlaneGeometry::guard_ring () {    return s_guard_ring;}
double   SiliconPlaneGeometry::ssd_gap () {       return s_ssd_gap;}
double   SiliconPlaneGeometry::ladder_gap () {    return s_ladder_gap;}
unsigned int SiliconPlaneGeometry::strips_per_die () {return s_stripPerWafer;}
unsigned int SiliconPlaneGeometry::n_si_dies() {return s_n_si_dies;}



/// insideActiveArea - does a point fall on (in local coords) an 
/// active part of the detector? (-) if no
/// spacing is the ladder gap for local X; the ssd gap for local Y
double SiliconPlaneGeometry::insideActiveArea1D (double x, double spacing) {
    x += panel_width()/2.;
    
    if (x < 0)
        return x;
    
    double edge = fmod(x, (die_width() + spacing));
    
    if (edge > die_width()) return -edge;  // edge is in the gap btwn si dies
    
    if (edge >= 0) {
        if (edge > die_width()/2.) edge -=die_width();
        edge = fabs(edge);
        edge -= guard_ring();
    }
    return edge;
}

/// insideActiveArea - does a point fall on (in local coords) an active
/// part of the detector? (-) if no
double SiliconPlaneGeometry::insideActiveArea (const HepPoint3D& p) {
    // Point q(p);
    // q.transform(globalToLocal());
    return insideActiveArea(p.x(), p.y());
}

/// stripId - compute the strip index corresponding to a given point
unsigned int SiliconPlaneGeometry::stripId (double x) {
    if (insideActiveArea(x, 1.0) > 0.) {    // check that this is inside the die
        x += panel_width()/2.;
        double  die = floor(x/(die_width() + ladder_gap())); // find die #
        double  wafer_x = fmod(x,(die_width() + ladder_gap())) - guard_ring();
        //if ((x/die_width() - die) >= 0.5)    die += 1.;
        
        return static_cast<unsigned>(floor(wafer_x/si_strip_pitch()) + 
            die*strips_per_die());
    } else return 65535;
}


/// localX - convert a strip ID to local coordinates
double SiliconPlaneGeometry::localX (unsigned ix) {
    double  nstrips = static_cast<double>(strips_per_die());
    double  i = fmod(static_cast<double>(ix),nstrips);
    // strip relative to wafer
    double  n = (die_width() + ladder_gap()) * 
        floor(static_cast<double>(ix)/nstrips);
    // die number
    
    return  n + (i + 0.5) * si_strip_pitch() - panel_width()/2. + guard_ring();
}


/// n_si_strips - return the number of Si strips in a single layer
unsigned int SiliconPlaneGeometry::n_si_strips ()  {
    return (n_si_dies() * strips_per_die());
}

/// si_strip_pitch - return the width of a single silicon strip
double SiliconPlaneGeometry::si_strip_pitch () {
    return ((die_width() - 2. * guard_ring()) / strips_per_die()); 
}

/// panel_width - width of all connected silicon in a single layer
double SiliconPlaneGeometry::panel_width () {
    return (n_si_dies() * die_width() + (n_si_dies() - 1) * ladder_gap());
}
