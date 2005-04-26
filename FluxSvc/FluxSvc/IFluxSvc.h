/** 
* @file IFluxSvc.h
* @brief definition of the interface for IFluxSvc
*
*  $Header$
*/
#ifndef _H_IFluxSvc
#define _H_IFluxSvc

// includes
#include "GaudiKernel/IInterface.h"
#include "CLHEP/Vector/Rotation.h"
#include "facilities/Observer.h"

#include <string>
#include <list>
#include <vector>
#include "astro/GPS.h"

// Declaration of the interface ID ( interface id, major version, minor version) 
static const InterfaceID IID_IFluxSvc("FluxSvc", 5 , 0); 

// forward declarations
class IFlux;
class IParticlePropertySvc;
class ISpectrumFactory;
class HepRandomEngine;

/** 
* \class IFluxSvc
* \brief The FluxSvc gaudi service interface
*
* \author Toby Burnett tburnett@u.washington.edu
* 
* $Header$
*/
class  IFluxSvc : virtual public IInterface {
public:

    /// just get an IFlux object by name -- see compositeSource
    virtual  StatusCode source(std::string name, IFlux*&)=0;


    /// return a list of legal names
    virtual std::list<std::string> fluxNames()const=0;

    /// add a new source
    virtual void addFactory(std::string name, const ISpectrumFactory* factory )=0;

    /// return pointer to the random engine that FluxSvc uses
    virtual HepRandomEngine* getRandomEngine()=0;

    /// pass a specific amount of time
    virtual void pass (double t)=0;    

    /// create a set of display windows using rootplot.
    virtual void rootDisplay(std::vector<const char*> arguments)=0;

    /// allow an external observer to attach to gps
    virtual void attachGpsObserver(Observer* anObserver)=0;

    /// Retrieve interface ID
    static const InterfaceID& interfaceID() { return IID_IFluxSvc; }

    ///return the pointer to the current IFlux object
    virtual IFlux* currentFlux()=0;

    /// name of the flux
    virtual std::string fluxName()const=0;

    /// set the glast tilt angles - this means glast will ALWAYS rock to these angles.
    virtual void setExplicitRockingAngles(double ang1,double ang2)=0;

    /// get the angular values of the satellite
    virtual std::pair<double,double> getExplicitRockingAngles()=0;

    /// set the desired pointing history file to use:
    virtual void setPointingHistoryFile(std::string fileName)=0;

    ///this transforms glast-local (cartesian) vectors into galactic (cartesian) vectors
    virtual HepRotation transformGlastToGalactic(double time)const=0;

    virtual HepRotation transformToGlast(double seconds,GPS::CoordSystem index)const=0;
    /// get the current satellite location
    virtual std::pair<double,double> location()=0;

    /// return pointer to the GPS instance of FluxSVc
    virtual GPS* GPSinstance()=0;

    ///set the type of rocking, along with the angle to rock to.
    ///0 = NONE : No rocking rotation done at all.
    ///1 = UPDOWN : Satellite will be rocked toward the north pole in the northern hemisphere, opposite in the south.
    ///2 = SLEWING : (experimental) like UPDOWN, except that rotation at equator happens gradually.
    ///3 = ONEPERORBIT : (needs work) LAT rocked northward for one orbit, southward for the next.
    ///4 = EXPLICIT :  Explicit angles given - this should be used with the setExplicit... function.
    ///5 = POINT:  Explicit pointing direction given - setExplicitRockingAngles are (l,b).
    ///6 = HISTORY - Filename given to stand for a pre-recorded pointing history.  Use the setPointingHistoryFile function.
    virtual std::vector<double> setRockType(int rockType, double rockAngle)=0;

    ///this should return the source file names, along with the contained sources.
    virtual std::vector<std::pair< std::string ,std::list<std::string> > > sourceOriginList() const=0;

    /// set pointer to a flux object, constructed from set of names
    virtual StatusCode compositeSource(std::vector<std::string> names, IFlux*& flux)=0;

    /// return current SAA status
    virtual bool insideSAA()=0;

};

#endif  // _H_IFluxSvc
