/** 
* @file SurfaceMuons.cxx
* @brief declaration and definition of SurfaceMuons
*
*  $Header$
*/
#include "Spectrum.h"
#include "SpectrumFactory.h"

#include "CLHEP/Random/RandFlat.h"
#include <string>
#include <utility>
#include <algorithm>
#include <map>

/** 
* \class SurfaceMuons
*
* \brief Spectrum representing cosmic ray muon flux at the Earth's surface
* \author T. Burnett
* 
* $Header$
*/
//

class SurfaceMuons : public Spectrum
{
public:
    /** @brief ctot
    @param paramstring string from xml
    */
    SurfaceMuons(const std::string& paramstring);

    /// flux integrated over energy, will be multipled by solidAngle to get a rate/m^2
    double flux (double time ) const { return m_flux;}

    double solidAngle()const{return 1.0;}
    
    /** @brief sample a single particle energy from the spectrum
    @param current time (ignored)
    @return the energy in GeV, or MeV if the source was specified as such
    */
    double energy( double time);

    
    /** @brief return solid angle pair (costh, phi) for the given energy
    @param the energy previously generated
    */
    virtual std::pair<double,double> dir(double energy);
    
    
    virtual std::string title() const{return "SuraceMuons";}
    virtual const char * particleName() const;
    inline  const char * nameOf() const {return "SurfaceMuons";}
    
    
private:

    // local function that approximates the spectrum as a function of E*cos(theta)
    double spectrum(double e);

    /**
    This is the integral spectrum as a map where the key is the integral spectrum and the
    value the associated energy. This makes is easy to invert it by using the STL algorithm lower_bound.

  */
    std::map<double,double> m_ispec; // integral spectrum, for sampling
    double m_flux;
    double m_index;
    double m_E0, m_emax;
    double m_costh;
    double m_total;
    bool m_vertical;
    
};


static SpectrumFactory<SurfaceMuons> factory;
const ISpectrumFactory& SurfaceMuonsFactory = factory;

SurfaceMuons::SurfaceMuons(const std::string& paramstring)
: m_index(2.71)
, m_flux(70.*2*M_PI/3) // 70 vertical, cos2theta after that
, m_emax(1000)
, m_E0(5) 
, m_total(0)
, m_vertical(false)
{
    //Purpose: Initializes parameters during construction
    
    std::vector<float> params;
    
    parseParamList(paramstring,params);


    m_vertical = (params.size()>0) && params[0]==-1 ; 

    // create integral table of the flux function, as a map of
    // energy and e*flux(e), with logarithmic energies
    int n=100;
    double  emin=1;
        
    for( int i=0; i< n; ++i){
        double 
            e = pow(10, 0.025*i),
            f= e*spectrum(e);
        
        m_ispec[m_total +=f] =e;
    }
}


const char* SurfaceMuons::particleName()const
{
/// purpose: return a point to the particle name, either mu+ or mu-
    static const char * pnames[] = {"mu+", "mu-"};
    static double 
        charge_ratio = 1.2,  // from many measurements.
        plus_fraction=1/(1+charge_ratio);
    return RandFlat::shoot()>plus_fraction? pnames[0]:pnames[1];
}
  
/// sample a single particle energy from the spectrum: assume called first
double SurfaceMuons::energy( double )
{
    using namespace std;
    // first choose the angle for the dir function
    // if doing vertical, don't!
    static double third=1.0/3.0;
    m_costh = m_vertical? 1 : pow(RandFlat::shoot(), third);

    // select an energy*costh from the distribution

    double fraction = RandFlat::shoot();
     map<double,double>::const_iterator 
         element  = m_ispec.lower_bound(fraction*m_total),
         prev = element;
         prev--;
     if( element!=m_ispec.end() ){
         double energy = prev->second;
         // interpolate if not last one
         double 
             deltaI = fraction-prev->first/m_total,
             deltaE = element->second -energy,
             correct = deltaE*deltaI;
         return (energy+deltaI)/m_costh; 
    }
    return m_emax/m_costh;
}

  
std::pair<double,double> SurfaceMuons::dir(double)
{
    // purpose: return the pair; note uses the previous value for cos(theta)
    return std::make_pair(m_costh, RandFlat::shoot(2*M_PI));
}

static inline double sqr(double x){return x*x;}

double SurfaceMuons::spectrum(double ecth)
{
    // purpose: spectrum as a function of E*cos(theta)
    // returns value differential in E, according to stuff in particle properties

    double 
        atmos = 1/(1+1.1*ecth/115.) + 0.054/(1+1.1*ecth/850.),
        cutoff = ecth<30? exp(-sqr(::log(ecth/30.))/0.55) : 1.0;


    return pow(ecth, -2.71)*atmos*cutoff;
}
