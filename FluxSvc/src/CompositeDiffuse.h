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

#include "dom/DOM_Document.hpp"
#include "dom/DOM_Element.hpp"
#include "xml/Dom.h"
#include "xml/IFile.h"
#include "xml/XmlParser.h"
#include <string>
#include <typeinfo>

//#include <cmath>
//#include <algorithm>
//#include <functional>
#include <fstream>
//#include <iostream>

class FluxSource;

//! holds multiple Eventsource objects ; acts as a container for them.
class CompositeDiffuse: public CompositeSource { 
public:
    
    typedef struct{
        double l;
        double b;
        double flux;
        double energyIndex;
    }PointSourceData;

    CompositeDiffuse(double totalFlux):
      m_totalFlux(totalFlux),m_unclaimedFlux(totalFlux)
      {
          fillTable();
          setFileFlux();
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
      
      ///fills the log N/log S characteristic table with data from the XML file
      void fillTable();
      
      //stolen from FluxMgr, this collects xml and dtd and prepares an input
      //file for parsing by the DOM parser.
      std::string writeXmlFile(const std::vector<std::string>& fileList);


      /// write the characteristics of the current source distribution to a stream
      void writeSourceCharacteristic(std::ostream& out){
        out << fullTitle() << std::endl;
        out << "Diffuse Added Point Sources:" << std::endl;
        out << "l" << '\t' << "b" <<'\t' << "flux" << '\t'<< "energyIndex" << std::endl << std::endl;

        std::vector<PointSourceData>::iterator now= m_listOfDiffuseSources.begin();
        for ( ; now != m_listOfDiffuseSources.end(); now++) {
            out << (*now).l << '\t' <<(*now).b << '\t' << (*now).flux << '\t' << (*now).energyIndex << std::endl;
        }

    }
      
private:
    double m_totalFlux; // The total flux from the entire sky, particles/sec/m^2/steradian.
    double m_unclaimedFlux; // The amount of the flux "unaccounted for" by the known sources
    double getRandomFlux();
    //long double pofi(long double intensity);
    long double logNlogS(long double flux);
    std::vector<std::pair<double,double> > m_logNLogS; // inputted logN.logS graph, to be used for finding random fluxes.
    std::vector<PointSourceData> m_listOfDiffuseSources; //information on the sources which got added randomly.
    double m_minFlux; //the minimum flux to be considered in the spectrum
    double m_maxFlux; //the maximum flux to be considered.
    double m_totalIntegratedFlux; //integrated flux, using input file.
    std::string m_dtd; //the DTD document for the XML data to be read in.
    /// list of sources for easy lookup
    std::map<std::string, DOM_Element > m_sources;
};
