
#include "MapSpectrum.h"

#include <cmath>
#include <algorithm>
#include <functional>
#include <fstream>
#include <iostream>

#include "fitsio.h"  
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
    /*
    //now get the second element from params to be the filename
    std::string input = params.c_str();
    int i=params.find_first_of(",");
    input = input.substr(i+1);  
    
    m_particle_name = "gamma";
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
    }*/
 

    //----------------------fitsio----------------------------------------------------
    int status=0;
    const char* flux_root = ::getenv("FLUXSVCROOT");
    //define the file
    std::string doc_path2= (flux_root? std::string(flux_root)+"/" : "");
    std::string fileName2 = doc_path2+"doc/bremss_skymap_41p_222111";
    //open the file
    fitsfile *fptr;
    ffopen(&fptr, fileName2.c_str(), 0, &status);
    printf("\Opening file: %s status = %d\n", fileName2.c_str(),status);

    int nfound, anynull;
    long naxes[2], fpixel, nbuffer, npixels, ii;
    long jj=0;

#define buffsize 720
    float datamin, datamax, nullval, buffer[buffsize];
    status = 0;

    /* read the NAXIS1 and NAXIS2 keyword to get image size */
    if ( fits_read_keys_lng(fptr, "NAXIS", 1, 2, naxes, &nfound, &status) )
        std::cout << "error: " << status << std::endl;

    npixels  = naxes[0] * naxes[1];         /* number of pixels in the image */
    fpixel   = 1;
    nullval  = 0;                /* don't check for null values in the image */
    datamin  = 1.0E30;
    datamax  = -1.0E30;

    while (npixels > 0)
    {
                double curl,curb,curint,curind;
      nbuffer = npixels;
      if (npixels > buffsize)
        nbuffer = buffsize;     /* read as many pixels as will fit in buffer */

      /* Note that even though the FITS images contains unsigned integer */
      /* pixel values (or more accurately, signed integer pixels with    */
      /* a bias of 32768),  this routine is reading the values into a    */
      /* float array.   Cfitsio automatically performs the datatype      */
      /* conversion in cases like this.                                  */

      if ( fits_read_img(fptr, TFLOAT, fpixel, nbuffer, &nullval,
                  buffer, &anynull, &status) )
        std::cout << "error: " << status << std::endl;

      for (ii = 0; ii < nbuffer; ii++)  {
          //std::cout << buffer[ii] << "  ";
          curb=jj/2-90;
          curl=ii/2-180;
          curint=buffer[ii];
          curind = 2.;
          m_catalog[std::make_pair<int,int>((int)curl,(int)curb)].intensity = curint;
          m_catalog[std::make_pair<int,int>((int)curl,(int)curb)].index = curind;
//          std::cout << (int)curl << ", " << (int)curb << ", " << m_catalog[std::make_pair<int,int>(curl,curb)].intensity*100000. << std::endl;

      }
            jj++;
      //std::cout << std::endl;
      npixels -= nbuffer;    /* increment remaining number of pixels */
      fpixel  += nbuffer;    /* next pixel to be read in image */
    }

//            input_file >> curl;
//            input_file >> curb;
//            input_file >> curint;
//            input_file >> curind;
//            m_catalog[std::make_pair<int,int>(curl,curb)].intensity = curint;
//            m_catalog[std::make_pair<int,int>(curl,curb)].index = curind;

    setNetFlux();
    std::cout << "net flux is " << m_netFlux << std::endl;
    if ( fits_close_file(fptr, &status) )
        std::cout << "error: " << status << std::endl;

//////////////////////////////////////////////////////////////////////////////

//    ffclos(fptr, &status);
//std::cout << "here" << std::endl;
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
    //double sizeof1by1 = 1.;
    double l=-180;
    double b=-90;
    while(remainingFlux>=0){
        //std::cout << "remainingFlux = " << remainingFlux;
        remainingFlux -= m_catalog[std::make_pair<int,int>(l,b)].intensity*sizeOf1by1(b);
        b++;
        if(b==90){b=-90; l++;}
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
    for(int l=-180 ; l<=180 ; l++){
        //double sizeof1by1 = 1.;
        for(int b=-90 ; b<=90 ; b++){
            m_netFlux+=m_catalog[std::make_pair<int,int>(l,b)].intensity*sizeOf1by1(b);
            //std::cout << (int)l << ", " << (int)b << ", " << m_catalog[std::make_pair<int,int>(l,b)].intensity*10000000. << std::endl;
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

double MapSpectrum::sizeOf1by1(double b){
return cos(b*M_PI/180.);
}