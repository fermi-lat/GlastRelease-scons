#ifndef MapSpectrum_H
#define MapSpectrum_H
// $Heading:$

#include "Spectrum.h"
#include "facilities/Observer.h"
#include <vector>
#include <string>
#include <map>

typedef struct{
    double intensity;
    double index;
}cellinfo;


//!  Spectrum that reads its differential spectrum from a table
class MapSpectrum : public Spectrum
{
public:
    /// params is the filename to read
    MapSpectrum(const std::string& params);
    
    
    /// return total flux 
    virtual double flux() const;
    
    
    /// sample a single particle energy from the spectrum
    virtual float operator() (float)const;
    
    ///this returns a galactic direction, in the form (l,b)
    std::pair<double,double> dir(double energy, HepRandomEngine* engine);
    
    
    virtual std::string title() const;
    virtual const char * particleName() const;
    inline  const char * nameOf() const {return "MapSpectrum";}
    //   use default destructor, copy constructor, and assignment op.
    
    double sizeOf1by1(double b);

    ///sets the net flux
    void setNetFlux();   
    float parseParamList(std::string input, int index);
    
private:
    std::map<std::pair<int,int>,cellinfo> m_catalog;
    
    double m_flux;   // current flux (set when direction changes)
    
    double m_index;  //current power-law index for the energy (set when direction changes)
    double m_netFlux; //sum over the entire sky's flux
    float m_E0;  //energy base
    
    
    typedef std::pair<double,double> efpair;
    std::vector<efpair> integ_flux;
    
    
    std::string initialization_document;
    std::string m_particle_name;
};



#endif // MapSpectrum_H