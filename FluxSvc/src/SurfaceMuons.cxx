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
    @param number between 0 and 1, probably generated randomely to sample spectum
    @return the energy in GeV, or MeV if the source was specified as such
    */

    virtual float operator() (float r)const;

    
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

    std::map<double,double> m_ispec; // integral spectrum, for sampling
    double m_flux;
    double m_index;
    double m_E0, m_emax;
    mutable double m_costh;
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
    int n=120;
    double  emin=1;
        
    for( int i=0; i< n; ++i){
        double 
            e = pow(10, 0.025*i),
            f= e*spectrum(e);
        
        m_ispec[e]= (m_total +=f);
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
float SurfaceMuons::operator() (float r)const
{
    // first choose the angle for the dir function
    // if doing vertical, don't!
    static double third=1.0/3.0;
    m_costh = m_vertical? 1 : pow(r, third); // const is a pain

    // select an energy*costh from the distribution

    for( std::map<double,double>::const_iterator it = m_ispec.begin(); 
    it!= m_ispec.end(); ++it) {

        double part=it->second;
        if(r*m_total > it->second) continue;
            
        // reached the level: return corresponding energy (TODO: adjust)
        return it->first/m_costh; 
    }
    return m_emax;
}

  
    /// return solid angle pair (costh, phi) 
std::pair<double,double> SurfaceMuons::dir(double)
{
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
