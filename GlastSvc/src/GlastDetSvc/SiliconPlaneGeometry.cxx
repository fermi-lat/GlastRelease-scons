// File and Version Information:
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

// static variable implementations, initialized from geometry
unsigned int SiliconPlaneGeometry::s_n_si_dies = 0;        // was 4; // # of Si dies across single SSD plane
unsigned int SiliconPlaneGeometry::s_stripPerWafer = 0;    // was 384;
double       SiliconPlaneGeometry::s_dice_width     = 0.0; // was 89.500; // width of a single silicon die 

// width of the dead area around the perimeter of a die (from edge to active)
double       SiliconPlaneGeometry::s_guard_ring     = 0.0; // was 0.974;

// width of the gap between SSD dies (glue & wirebonds, etc.) (cm)
double       SiliconPlaneGeometry::s_ssd_gap        = 0.0; // was 0.025;
double       SiliconPlaneGeometry::s_ladder_gap     = 0.0; // was 0.200;


double SiliconPlaneGeometry::insideActiveArea (double x, double y) 
{
    // Purpose and Method: does a point fall on (in local coords) an active
    //   part of the detector? (-) if no
    
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

double SiliconPlaneGeometry::insideActiveArea1D (double x, double spacing) {
    // Purpose and Method: does a point fall on (in local coords) an
    //  active part of the detector? (-) if no
    //  spacing is the ladder gap for local X; the ssd gap for local Y
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

double SiliconPlaneGeometry::insideActiveArea (const HepPoint3D& p) {
    // Purpose and Method:  does a point fall on (in local coords) an active
    //   part of the detector? (-) if no
    
    // Point q(p);
    // q.transform(globalToLocal());
    return insideActiveArea(p.x(), p.y());
}

unsigned int SiliconPlaneGeometry::stripId (double x) {
    // Purpose and Method:  compute the strip index corresponding to a 
    //  given point
    
    if (insideActiveArea(x, 1.0) > 0.) {    // check that this is inside the die
        x += panel_width()/2.;
        double  die = floor(x/(die_width() + ladder_gap())); // find die #
        double  wafer_x = fmod(x,(die_width() + ladder_gap())) - guard_ring();
        //if ((x/die_width() - die) >= 0.5)    die += 1.;
        
        return static_cast<unsigned>(floor(wafer_x/si_strip_pitch()) + 
            die*strips_per_die());
    } else return 65535;
}

double SiliconPlaneGeometry::localXDouble (double x) {
    // Purpose and Method:  convert a non-integer strip ID to local coordinates
    double  nstrips = strips_per_die();
    double  i = fmod(x,nstrips);
    // strip relative to wafer
    double  n = (die_width() + ladder_gap()) * 
        floor(x/nstrips);
    // die number
    
    return  n + (i + 0.5) * si_strip_pitch() - panel_width()/2. + guard_ring();
}

double SiliconPlaneGeometry::localX (unsigned ix) {
    // Purpose and Method:  convert a strip ID to local coordinates
    double x = ix;
    return localXDouble(x);
    /*
    double  nstrips = static_cast<double>(strips_per_die());
    double  i = fmod(static_cast<double>(ix),nstrips);
    // strip relative to wafer
    double  n = (die_width() + ladder_gap()) * 
    floor(static_cast<double>(ix)/nstrips);
    // die number
    
      return  n + (i + 0.5) * si_strip_pitch() - panel_width()/2. + guard_ring();
    */
}

HepPoint3D SiliconPlaneGeometry::siPlaneCoord( const HepPoint3D& p, 
                                              idents::VolumeIdentifier id) {
    // Purpose and Method:  return a plane coordinate given a ladder coordinate
    
    double xShift = 0, yShift = 0;
    int view      = id[5];
    if (id.size()==9) {
        int ladderNum = id[7];
        //std::cout << " sgp ladder " << ladderNum << std::endl;
        xShift = 0.5*die_width() + ladderNum*(die_width() + ladder_gap()) - 0.5*panel_width();
        if (view==1) xShift = -xShift;
        
        double panelLength = n_si_dies()*die_width() + (n_si_dies()-1)*ssd_gap();
        yShift = 0.5*die_width() + ladderNum*(die_width() + ssd_gap()) - 0.5*panelLength;
    }
    
    xShift+= p.x();
    //std::cout << " spg view " << view << " " << xShift ;
    if (view==1) xShift = -xShift; 
    //std::cout << " " << xShift  << std::endl;
    yShift+= p.y();
    
    return HepPoint3D(xShift, yShift, p.z());
}


unsigned int SiliconPlaneGeometry::n_si_strips ()  {
    // Purpose and Method:  return the number of Si strips in a single layer
    return (n_si_dies() * strips_per_die());
}

double SiliconPlaneGeometry::si_strip_pitch () {
    // Purpose and Method:  return the width of a single silicon strip
    return ((die_width() - 2. * guard_ring()) / strips_per_die()); 
}

double SiliconPlaneGeometry::panel_width () {
    // Purpose and Method:  width of all connected silicon in a single layer
    return (n_si_dies() * die_width() + (n_si_dies() - 1) * ladder_gap());
}

void SiliconPlaneGeometry::init(IGlastDetSvc *p_detSvc) 
{	
    StatusCode sc;
    double temp;
    sc = p_detSvc->getNumericConstByName("ssdGap", &s_ssd_gap);
    sc = p_detSvc->getNumericConstByName("ladderGap", &s_ladder_gap);	
    sc = p_detSvc->getNumericConstByName("stripPerWafer", &temp);
    s_stripPerWafer = temp + 0.5;
    sc = p_detSvc->getNumericConstByName("nWaferAcross", &temp);
    s_n_si_dies = temp + 0.5;
    sc = p_detSvc->getNumericConstByName("SiWaferSide", &s_dice_width);	
    double SiActiveSide;
    sc = p_detSvc->getNumericConstByName("SiWaferActiveSide", &SiActiveSide);	
    s_guard_ring = 0.5*(s_dice_width - SiActiveSide);
    
    /*std::cout << "SPG init: ssdgap " << s_ssd_gap << " laddergap " << s_ladder_gap 
    << " stripperwafer " << s_stripPerWafer << " ndies " << s_n_si_dies 
    << " waferside " << s_dice_width << " deaddist " << s_guard_ring << std::endl;
    */
    
    return;
}
