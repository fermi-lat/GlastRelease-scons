//	$Header$

#ifndef FluxSource_h
#define FluxSource_h 1
/** 
* \class FluxSource
*
* \brief EventSource subclass to take over the functionality of the old Flux class, 
* which implemented a GISMO based event generation scheme.
* 
* $Header$
*/

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
    FluxSource ( ISpectrum* aSpec = 0, double aFlux = 0 );
    FluxSource ( const DOM_Element& xelem );
    FluxSource::FluxSource(double aFlux, ISpectrum* aSpec,  double l, double b);
    
    ///    destructor
    virtual ~FluxSource();
    
    ///    generate an event from a Flux object ??
    virtual FluxSource* event(double time);
    
    ///    full-length title description of this EventSource.
    virtual std::string fullTitle () const;
    
    ///    brief title description (for display) for this event source
    virtual std::string displayTitle () const;
    
    ///    getLaunch - compute launch point, direction, & energy
    virtual void computeLaunch (double time=0);
    
    virtual double flux(double time)const; // calculate flux for attached spectrum
    
    /// return effective solid angle
    double solidAngle()const;

    /// return a title describing the spectrum and angles
    std::string title()const;
    
    /// print facility
    void  printOn ( std::ostream&  ) {}
    
    /// set spectrum, with optional parameter to set the maximum energy?
    void spectrum(ISpectrum* s, double emax=-1);
    
    ISpectrum* spectrum() const{ return m_spectrum; }

    double interval(double ){return m_interval;}

    //! Denotes what Energy Units the energy
    //! of incoming particles are in
    enum EnergyScale { 
        MeV,        //! MeV
        GeV         //! GeV
    } m_energyscale;
    virtual int eventNumber()const;
    
    double energy()const { return m_energy;}
    const HepVector3D& launchDir()const {return m_correctedDir;}//m_correctForTilt*m_launchDir;}
    const HepPoint3D&  launchPoint()const { return m_launchPoint;}
    
    static	double	s_radius;

  private:
  
      // base class for direction and position strategy
      class LaunchStrategy;

      // forward declaration of classes that handle the lauch direction
      class LaunchDirection;  // base class
      class RandomDirection;  // choose randomly from range 
      class SourceDirection;  // choose from an external source class
   
      // forward declaration of classes that handle launch point
      class LaunchPoint;  // base class
      class RandomPoint; // random strategy
      class FixedPoint;  // fixed, or pencil
      class Patch;  // a box

 
      LaunchPoint* m_launch_pt; // pointer to actual point stategy: must be set
      LaunchDirection* m_launch_dir;
      
      ISpectrum*         m_spectrum;	    // spectrum to generate

      double m_energy;
      // associated with a specific launch

      ///use GPS to correct m_launchDir for the rocking of the spacecraft.
      void correctForTiltAngle();
      
      /// result of strategy
      HepVector3D m_launchDir;

      ///direction after being corrected for the "tilt" angles.
      HepVector3D m_correctedDir;
      
      HepPoint3D  m_launchPoint;

      //!the "extra time" a source needs to come out of occlusion.
      double m_extime;
      double m_interval; //the current value of the interval in time to the next particle.

      //! whether or not the current particle is occluded by the earth
      bool occluded();    
      double calculateInterval (double time);

      ///interval function to be used by non-spectrum sources
      double explicitInterval (double time);

//      double m_maxEnergy; // max kinetic energy allowed when running a spectrum

};
#endif
