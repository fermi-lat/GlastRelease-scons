#if !defined(AFX_ISPECTRUM_H__9D00078B_7121_473C_813A_827B2F4126A5__INCLUDED_)
#define AFX_ISPECTRUM_H__9D00078B_7121_473C_813A_827B2F4126A5__INCLUDED_

/** 
* \class ISpectrum
*
* \brief The virtual interface for Spectrum-type objects.
*
* Class for holding function definitions of Spectrums...
* an abstract base class
* \author Sean Robinson
*
* $Header$
*/

#include <string>
#include <utility> // for std::pair


class ISpectrum  
{
    
public:
    
    ///  particle name that must be known to the particle service
    virtual const char * particleName()const=0;
    
    /** calculate the flux, particles/m^2/sr. (default zero)
        @param time the mission elapsed time in sec
    */
    virtual double    flux (double time ) const=0;
    
    /// return effective solid angle that will be used to determine the actual rate 
    virtual double solidAngle()const{return  1.0;} //flag that doesen't calculate.
    
    /// return a title describing the spectrum	
    virtual std::string title()const=0;
    
    /*! a (randomized) interval to the next event.  
    @param time the mission elapsed time in sec
     For time-independent rate, should correspond to exponential( 1/rate() )
     Return negative to do this with flux()*solidAngle().
     needs to know the cross-sectional area?
    */
    virtual double interval (double time)=0;
    
    /// return energy, either GeV or MeV
    virtual double energy( double time=0)=0;

    /** return direction in a pair:
    @param energy The generated energy, from previous call to energy
    @return direction is either in the format (cos theta, phi) for
    (zenith-local coordinates, or (l,b) (galactic coordinates).
    */
    virtual std::pair<double,double> dir(double energy)=0;
    
};


#endif // !defined(AFX_ISPECTRUM_H__9D00078B_7121_473C_813A_827B2F4126A5__INCLUDED_)
