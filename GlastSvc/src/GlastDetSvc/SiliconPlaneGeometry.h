//$Header$

#ifndef SiliconPlaneGeometry_H
#define SiliconPlaneGeometry_H

#include "CLHEP/Geometry/Vector3D.h"
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
/**
    Static functions only for geometry of individual Si plane
  */
class SiliconPlaneGeometry {
public:    
    /// (-) if non-active or gap (P is in local system)    
    static double insideActiveArea (const HepPoint3D& p);  
      
    /// compute strip id from local coordinate
    static unsigned int stripId (double x);

   /// location of strip ix in local coords, allows for non-integer strip number
    static double localXDouble ( double x);    
     
    /// location of strip ix in local coords
    static double localX ( unsigned int ix); 
	  
    /// insideActiveArea for local coords
    static double insideActiveArea(double x, double y); 
    
    /// insideActiveArea for local coord (one-D)
	/// allows for different spacing in X and Y
    static double insideActiveArea1D(double x, double spacing);

    /// number of silicon dies across a single layer
    static unsigned int n_si_dies ();
    
    /// number of silicon strips across a single die
    static unsigned int strips_per_die ();
    
    /// width of a single si die
    static double die_width ();
    
    /// number of silicon strips in a single layer
    static unsigned int n_si_strips ();
    
    /// number of front-end chips on a single si die
    static unsigned int n_fe_chips ();
    
    /// silicon strip pitch (computed from n_fe_chips & active are of the die)
    static double si_strip_pitch ();
    
    
    /// gap between the edge of the silicon die and the tray closeout
    static double electronics_gap ();
    
    /// gap between Silicon dies
    static double ssd_gap ();
    
    /// gap between ladders
    static double ladder_gap ();
    
    /// width of dead area, aka guarg_ring
    static double guard_ring ();
    static double panel_width(); 

	static void init(IGlastDetSvc * p_detSvc);

private:

    static double s_dice_width;
    static double s_ssd_gap;
    static double s_ladder_gap;
    static double s_guard_ring;
    static double s_panel_width;
    static unsigned int s_stripPerWafer;
    static unsigned int s_n_si_dies;

	static double s_waferSide;
	static double s_waferActive;

};

#endif

