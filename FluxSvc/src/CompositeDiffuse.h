// CompositeDiffuse.h: interface for the CompositeDiffuse class.
/** 
* \class CompositeDiffuse
*
* \brief CompositeDiffuse functions like CompositeSource, with added point sources
* CompositeDiffuse takes not only a list of initial sources, but also a total flux 
* over the whole sky.  It then attempts to fill in the remaining flux by generating point sources according
* to a logN/logS characteristic.
* \author Sean Robinson, University of Washington, 2002
* 
* $Header $
*/

#include "FluxSvc/EventSource.h"
#include "CompositeSource.h"
#include <vector>

//#include <cmath>
//#include <algorithm>
//#include <functional>
#include <fstream>
//#include <iostream>

class FluxSource;

//! holds multiple Eventsource objects ; acts as a container for them.
class CompositeDiffuse: public CompositeSource { 
public:
    
    CompositeDiffuse(double totalFlux):
      m_totalFlux(totalFlux),m_unclaimedFlux(totalFlux)
      {
          double logIntensity; // log of source intensity, in units of particles/sec/m^2/steradian.
          double logSourcesPerFifthDecade; //log of number of sources at an intensity, integrated over 1/5 decade.
          //double intensity; //source intensity, in units of particles/sec/m^2/steradian.
          //double sourcesPerFifthDecade; //number of sources at an intensity, integrated over 1/5 decade.
          
          std::string initialization_document = "src/lognlogs.txt";
          int fitsflag = 1; //1 if we have a FITS file, 0 if not
          const char* flux_root = ::getenv("FLUXSVCROOT");
          //define the file
          std::string doc_path2= (flux_root? std::string(flux_root)+"/" : "");
          std::string fileName = doc_path2+initialization_document;//"doc/bremss_skymap_41p_222111";
          
          // get the input file name, and open it to read the logN/logS plot.
          std::ifstream input_file;
          input_file.open(fileName.c_str());
          
          if(false == input_file.is_open())
          {
              std::cerr << "ERROR:  Unable to open:  " << fileName.c_str() << std::endl;
              exit(0);
          }
          else
          {
              //double curl,curb,curint,curind;
              
              while (!input_file.eof()){
                  //get the info in the file and plug it into 
                  input_file >> logIntensity;
                  input_file >> logSourcesPerFifthDecade;
                  std::cout << logIntensity << "   " << logSourcesPerFifthDecade<<std::endl;
                  //intensity = pow(10,logIntensity);
                  //sourcesPerFifthDecade = pow(10,logSourcesPerFifthDecade);
                  m_logNLogS.push_back(std::make_pair<double,double>(logIntensity, logSourcesPerFifthDecade));
              }
              //now figure out the maximum and minimum fluxes.
              m_minFlux = pow(10,(*m_logNLogS.begin()).first);
              m_maxFlux = pow(10,(m_logNLogS.back()).first);
              //std::cout << (*m_logNLogS.begin()).first << "   " << (m_logNLogS.back()).first <<std::endl;
              setFileFlux();

          }
          }
          
          ~CompositeDiffuse(){}
          
          
          /// generate an event from from one of the sources 
          /// which make up the composite, and return a pointer to it
          FluxSource* event (double time);
          
          /// Same as from CompositeSource, but allowing for total flux adjustment.
          void CompositeDiffuse::addSource (EventSource* aSource);
          
          ///Randomly determines a new source, and adds it.
          void addNewSource();
          
          double remainingFluxInterval();
          
          //sets the total flux integrated over the input file, for later use
          void setFileFlux();
          
private:
    double m_totalFlux; // The total flux from the entire sky, particles/sec/m^2/steradian.
    double m_unclaimedFlux; // The amount of the flux "unaccounted for" by the known sources
    double getRandomFlux();
    //long double pofi(long double intensity);
    long double logNlogS(long double flux);
    std::vector<std::pair<double,double> > m_logNLogS; // inputted logN.logS graph, to be used for finding random fluxes.
    double m_minFlux; //the minimum flux to be considered in the spectrum
    double m_maxFlux; //the maximum flux to be considered.
    double m_totalIntegratedFlux; //integrated flux, using input file.
      };
