
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
    
    //now get the second element from params to be the filename
    std::string input = params.c_str();
    //output.push_back(f);
    int i=params.find_first_of(",");
    input = input.substr(i+1);  
    
    m_particle_name = "gamma";
    // TODO: have params be an index 
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
    std::cout << "opening file " << fileName.c_str() << " for reading" << std::endl;
    input_file.open(fileName.c_str());
    
    if(false == input_file.is_open())
    {
        std::cerr << "ERROR:  Unable to open:  " << fileName.c_str() << std::endl;
        exit(0);
    }
    else
    {
        //we have openend the file here.
//        const int buffer_size = /*256*/30;
//        char buffer[buffer_size];
        
        //get a line of the file and put it into the buffer
//        input_file.getline(buffer,buffer_size,'\n');

        double curl,curb,curint,curind;
        
        //FILE *fp;
        //float ptmp,p[20];
        //int i=0;
        //int iline=0; 
        while (!input_file.eof()/*input_file != '\n'*/){
            //scanf(buffer,"%f",&ptmp);
           input_file >> curl;
           input_file >> curb;
           input_file >> curint;
           input_file >> curind;
           m_catalog[std::make_pair<int,int>(curl,curb)].intensity = curint;
           m_catalog[std::make_pair<int,int>(curl,curb)].index = curind;
 
 //           p[i++]=ptmp;
 //           if (i==4){
 //               i=0;
 //               iline++;
 //               //hist2->Fill(p[0],p[1]);
 //               curl=p[0];
 //               curb=p[1];
 //               curint=p[2];
 //               curind=p[3];
 //               m_catalog[std::make_pair<int,int>(curl,curb)].intensity = curint;
 //               m_catalog[std::make_pair<int,int>(curl,curb)].index = curind;
 //           }
        }
        
        
        /*while(scanf(buffer,"%f",&curl)){
        //here we want to fill our table
        input_file >> curl;
        input_file >> curb;
        input_file >> curint;
        input_file >> curind;
        m_catalog[std::make_pair<int,int>(curl,curb)].intensity = curint;
        m_catalog[std::make_pair<int,int>(curl,curb)].index = curind;
    }*/
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
    while(remainingFlux>=0){
        while(remainingFlux>=0){
            remainingFlux -= m_catalog[std::make_pair<int,int>(l,b)].intensity*sizeof1by1;
            b++;
        }
        l++;
    }
    
    //now set the flux:
    m_flux = m_netFlux;
    
    //and the energy index of this particle:
    m_index = m_catalog[std::make_pair<int,int>(l,b)].index;
    
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
    //for(int l=-180 ; l<=180 ; l++){
    //    for(int b=-90 ; b<=90 ; b++){
    //        m_catalog[std::make_pair<int,int>(l,b)].intensity = 0;
    //    }
    //}
    for(int l=-180 ; l<=180 ; l++){
        double sizeof1by1 = 1.;
        for(int b=-90 ; b<=90 ; b++){
            m_netFlux+=m_catalog[std::make_pair<int,int>(l,b)].intensity*sizeof1by1;
        }
    }
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