// $Header$
// 
//  Original author: Toby Burnett tburnett@u.washington.edu

#ifndef _H_IFlux_
#define _H_IFlux_


// includes
#include <string>
#include "CLHEP/Geometry/Point3D.h"
#include "CLHEP/Geometry/Vector3D.h"
#include "../src/SpectrumFactoryTable.h"
#include "geometry/CoordTransform.h"

class ParticleProperty;

//!  Abstract interface for an object that generates particles, Flux
class IFlux {
public:
    /// ctor, select the name
    IFlux(std::string name=""){};
    virtual ~IFlux(){}
    
    /// name of the flux
    virtual std::string name()const=0;
    
    /// full title of the flux
    virtual std::string title()const = 0;
    
    /// generate a new entry trajectory
    virtual void generate()=0;
    
    /// the particle name of the last particle generated 
    virtual std::string particleName()const=0;

    /// the particle property entry for the last particle generated 
    //virtual ParticleProperty* property()const=0;
    
    /// its kinetic energy
    virtual double energy()const=0;
    
    /// starting point 
    virtual HepPoint3D launchPoint()const=0;
    
    /// direction
    virtual HepVector3D launchDir()const=0;
    
    /// time (s) (absolute or elapsed??)
    virtual double time()const=0;
    
    /// return rate ( /mm**2 /s)
    virtual double rate()const=0;
    
    /// set the area of the target
    virtual void setTargetArea( double area)=0;
    
    /// get the target area
    virtual double targetArea()const =0;
    
    /// find which spectrum created the current particle
    virtual std::string findSource()const=0;
    
    /// return a unique number correcponding to that spectrum
    virtual int numSource()const=0;

    /// pass a specific amount of time
    virtual void pass ( double t)=0;

    ///get the transformation matrix due to orientation of the Galaxy 
    virtual Rotation CELTransform(double time)const=0;

    ///get the transformation matrix due to orientation of the spacecraft.
    virtual Rotation orientTransform(double time)const=0;
       
    virtual void addFactory(std::string name, const ISpectrumFactory* factory )=0;

    virtual /*int*/double gpsTime()const=0;

#if 0
    // get a description of the parameters that can be modified, and reference to a list of them
    virtual std::string paramlist(std::vector<double>& params)=0; 
    
    // set the parameters
    virtual void setParams(std::vector<double>& params)=0;
#endif
};


#endif // _H_FluxSvc
