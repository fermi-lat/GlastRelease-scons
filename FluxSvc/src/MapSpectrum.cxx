
#include "MapSpectrum.h"

#include <cmath>
#include <algorithm>
#include <functional>
#include <fstream>
#include <iostream>

#include "SpectrumFactory.h"

static SpectrumFactory<MapSpectrum> factory;
const ISpectrumFactory& MapSpectrumFactory = factory;


MapSpectrum::MapSpectrum(const std::string& params)
:m_E0(parseParamList(params,0))
{
    //let's set the map to zero right off the bat:
    for(int l=-180 ; l<=180 ; l++){
        //double sizeof1by1 = 1.;
        for(int b=-90 ; b<=90 ; b++){
            m_catalog[std::make_pair<int,int>(l,b)].intensity = 0;
        }
    }
    //now get the second element from params to be the filename
    std::string input = params.c_str();
    //output.push_back(f);
    int i=params.find_first_of(",");
    input = input.substr(i+1);  
    
    m_particle_name = "gamma";
    //std::cout << "here" << std::endl;
    if(input.empty())
        initialization_document = "doc/test.txt";
    else
        initialization_document = input.c_str();
    
    // construct filename, assuming data files in folder <packageroot>/doc
    std::string fileName = "";
    const char* flux_root = ::getenv("FLUXSVCROOT");
    std::string doc_path= (flux_root? std::string(flux_root)+"/" : "");
    fileName = doc_path+initialization_document;
    std::ifstream input_file;
    //std::cout << "opening file " << fileName.c_str() << " for reading" << std::endl;
    input_file.open(fileName.c_str());
    
    if(false == input_file.is_open())
    {
        std::cerr << "ERROR:  Unable to open:  " << fileName.c_str() << std::endl;
        exit(0);
    }
    else
    {
        double curl,curb,curint,curind;
        
        while (!input_file.eof()){
            input_file >> curl;
            input_file >> curb;
            input_file >> curint;
            input_file >> curind;
            m_catalog[std::make_pair<int,int>(curl,curb)].intensity = curint;
            m_catalog[std::make_pair<int,int>(curl,curb)].index = curind;
        }
        setNetFlux();
    }
}


/// calculate flux for the current position
double MapSpectrum::flux() const
{
    return m_flux;
}

///this returns a galactic direction, in the form (l,b)
std::pair<double,double> MapSpectrum::dir(double energy, HepRandomEngine* engine){
    //here is where we decide where the next particle will come from, and set the associated flux and energy of the particle.
    double remainingFlux=RandFlat::shoot(m_netFlux);
    double sizeof1by1 = 1.;
    double l=-180;
    double b=-90;
    //while(remainingFlux>=0){
    while(remainingFlux>=0){
        //std::cout << remainingFlux << std::endl;
        remainingFlux -= m_catalog[std::make_pair<int,int>(l,b)].intensity*sizeof1by1;
        //std::cout << "l=" << l << ",b=" << b << std::endl;
        b++;
        if(b==90){b=-90; l++;}
    }
    //}
    
    //now set the flux:
    m_flux = m_netFlux;
    
    //and the energy index of this particle:
    m_index = m_catalog[std::make_pair<int,int>(l,b)].index;
    //std::cout << "l=" << l << ",b=" << b << std::endl;
    
    return std::make_pair<double,double>(l,b);    
}


/// sample a single particle energy from the spectrum
float MapSpectrum::operator() (float f)const
{
    double m_emax = 100.;
    if( m_index == 0.0 )     return m_E0;
    
    if( m_index == 1.0 ) return m_E0*exp(f*log(m_emax/m_E0));
    
    float x = 1 - exp((1-m_index)*log(m_emax/m_E0));
    return m_E0*exp(log(1-x*f)/(1-m_index));
}


std::string MapSpectrum::title() const
{
    return "MapSpectrum";
}

const char * MapSpectrum::particleName() const
{
    return m_particle_name.c_str();
}

void MapSpectrum::setNetFlux(){
    m_netFlux=0.;
    for(int l=-180 ; l<=180 ; l++){
        double sizeof1by1 = 1.;
        for(int b=-90 ; b<=90 ; b++){
            m_netFlux+=m_catalog[std::make_pair<int,int>(l,b)].intensity*sizeof1by1;
        }
    }
    //std::cout << m_netFlux << std::endl;
    return;
}

float MapSpectrum::parseParamList(std::string input, int index)
{
    std::vector<float> output;
    int i=0;
    for(;!input.empty() && i!=std::string::npos;){
        float f = ::atof( input.c_str() );
        output.push_back(f);
        i=input.find_first_of(",");
        input= input.substr(i+1);
    } 
    return output[index];
}