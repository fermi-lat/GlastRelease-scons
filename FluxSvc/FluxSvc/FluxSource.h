//	$Header$
//	EventSource subclass to take over the functionality of the old Flux class, which implemented
//	a GISMO based event generation scheme.

#ifndef FluxSource_h
#define FluxSource_h 1

#include "FluxSvc/EventSource.h"
#include "FluxSvc/ISpectrum.h"
#include "geometry/Point.h"

// forward declarations
class Box;
class DOM_Element;

// 
//! FluxSource:  class which manages to compute flux from various particle source configurations


class FluxSource : public EventSource  
{
  public:      
      ///  constructor
      FluxSource ( ISpectrum* aSpec = 0, double aFlux = /*1*/0 );
      FluxSource ( const DOM_Element& xelem );
      FluxSource::FluxSource(double aFlux, ISpectrum* aSpec,  double l, double b);

      ///    destructor
      virtual ~FluxSource();
          
      ///    generate an event from a Flux object ??//.
      virtual FluxSource* event(double time);
      
      ///    full-length title description of this EventSource.
      virtual std::string fullTitle () const;
      
      ///    brief title description (for display) for this event source
      virtual std::string displayTitle () const;
      
      ///    getLaunch - compute launch point, direction, & energy
      virtual void computeLaunch (double time=0);

      virtual double flux(double time)const; // calculate flux for attached specturm

      /// return effective solid angle
      double solidAngle()const;

      /// set parameters solid angle coverage and generation area
      void setPhiRange(double min_phi, double max_phi);
      void setCosThetaRange(double min_cos, double max_cos);

      /// set the mode for calculating direction to use the spectrum object. 
      /// it then must implement a solidAngle(double)const to return the effective
      /// solid angle 
      void useSpectrumDirection();

      void setAcceptance ();

      /// provide posibility of overriding box for specific launch direction and point.
      /// if only dir specified, will coose randomly
      void setLaunch(const Vector& dir, const Point& pos);
      void setLaunch(const Vector& dir);
      void setLaunch(double theta, double phi);
      void setLaunch(double xMax, double xMin, double yMax, double yMin, 
                     double zTop, double zBot, bool fan);
      void FluxSource::setLaunch(double theta, double phi, double xMax, 
                                 double xMin, double yMax, double yMin, 
                                 double zTop, double zBot);

      void unSetLaunch();
      
      // return of angle & solid angle parameters
      double minPhi() const      {return _minPhi;}
      double maxPhi() const      {return _maxPhi;}
      double minCosTheta() const {return _minCos;}
      double maxCosTheta() const {return _maxCos;}
      double phi() const         {return _phi;}
      double theta() const       {return _theta;}

      
      void getSurfacePosDir(); // setup position & direction for Surface illumination

      /// return a title describing the spectrum and angles
      std::string title()const;

      /// print facility
      void  printOn ( std::ostream& out ) {}
      
      /// set spectrum, with optional parameter to set the maximum energy?
      void spectrum(ISpectrum* s, double emax=-1);
      ISpectrum* spectrum() const{ return m_spectrum; }


      double calculateInterval (double time);
      double interval(double time){return m_interval;}

      void FluxSource::getGalacticDir(double l,double b);

      /// acess to the maximum energy (kinetic)
      double maxEnergy()const { return m_maxEnergy;}
      void setMaxEnergy(double e) { m_maxEnergy = e; }


      ///use GPS to correct m_launchDir for the rocking of the spacecraft.
      void FluxSource::correctForTiltAngle();

     //! choices for generating incoming particle trajectory
      enum LaunchType { 
          NONE,         //! random direction
          POINT,        //! fixed pointand direction(unused)
          DIRECTION ,   //! fixed direction
          SURFACE,      //! random point and direction within a fixed surface
          SPECTRUM, //! direction calculated by the spectrum object
          SPECGAL,   //! direction coming from the spectrum object, in the form (l,b) - galactic coordinates
          PATCHFIXED,   //! fixed direction, fixed surface (unused)
          GALACTIC  //!  fixed direction with respect to the galactic coordinate system (l,b)
      } m_launch;

	  enum PointType { 
          NOPOINT,         //! random point
          SINGLE,        //! fixed point
          PATCH   //! , fixed surface
      } m_pointtype;

          //! Denotes what Energy Units the energy
          //! of incoming particles are in
	  enum EnergyScale { 
          MeV,        //! MeV
          GeV         //! GeV
      } m_energyscale;

/*
          //! Denotes what coordinate system the coordinates
          //! of incoming particles are in
	  enum FrameType { 
          EARTH,        //! relative to the earth (default)
          GLAST,        //! relative to the glast satellite(differs from EARTH by the angle GLAST is rotated off zenith)
          GALAXY      //! relative to the galactic origin
      } m_frametype;
*/

 
          //! whether or not the current particle is occluded by the earth
          bool occluded();

      virtual int eventNumber()const;

      double energy()const { return m_energy;}
      const Vector& rawDir()const {return m_launchDir;}
      const Vector& launchDir()const {return m_correctedDir;}//m_correctForTilt*m_launchDir;}
      const Point&  launchPoint()const { return m_launchPoint;}

	  void refLaunch(LaunchType launch);
	  void refPoint(PointType point);
	  int refLaunch()const {return m_launch;}
	  int refPoint()const {return m_pointtype;}

  private:

      static	double	s_radius, s_backoff;

      ISpectrum*         m_spectrum;	    // spectrum to generate

      double m_maxEnergy; // max kinetic energy allowed when running a spectrum
      
      double _minCos, _maxCos, _minPhi, _maxPhi;
      // solid angle, with respect to z-axis
      
      double m_rmin, m_rmax;
      // minimum, maxumum value for random number (fraction of spectrum)
      
      double  _phi, _theta;
      //  angles, with respect to z-axis
      
      // galactic position of source - for use with GALACTIC
      double m_galb, m_gall;
      // dimensions of patch to be illuminated  -- for use with SURFACE
      double patchHeight, patchBottom, patchTop, patchXmax, patchXmin;
      double patchYmax, patchYmin, patchWidX, patchWidY, Fratio;
      double patchRange, patchOffset;
      // for use with SURFACE -- sidePatch means illumination of ONE side
      // fanBeam == false means polar beam
      bool sidePatch, fanBeam;
      Box* illumBox;
      double m_interval; //the current value of the interval in time to the next particle.
      

      //! transform the current m_launchDir into GLAST-relative coordinates
      //! to give to the client
      //void transformDirection();

      Vector m_launchDir;
      
      //Vector m_transformDir;
      Point  m_launchPoint;
      double m_energy;
      // associated with a specific launch

      /// rotation associated with the "tilting" angles.
      Rotation m_correctForTilt;

      ///direction after being corrected for the "tilt" angles.
      Vector m_correctedDir;
      //!the "extra time" a source needs to come out of occlusion.
      double m_extime;

    void randomLaunchPoint(); 
    // calculate a random launch point in a plane perpendicular to _launchDir

    ///interval function to be used by non-spectrum sources
    double FluxSource::explicitInterval (double time);
};
#endif
