// CompositeDiffuse.h: interface for the CompositeDiffuse class.
/** 
* \class CompositeDiffuse
*
* \brief CompositeDiffuse functions like CompositeSource, with added point sources
* CompositeDiffuse takes not only a list of initial sources, but also a total flux 
* over the whole sky.  It then attempts to fill in the remaining flux by generating point sources according
* to a logN/logS characteristic.
* First:  the constructor calls the functions to read in the logN/logS characteristic
* from the xml datatable.  Note:  this characteristic is assumed to be in non-differential format
* (i.e. N is the number of sources with flux equal to or higher than S).
* Second:  the Differential logN/logS curve is obtained and integrated over to find the total unresolved flux.
* m_currentremainingFlux holds (flux, sources in bin) information.
* Third:  each time a particle is needed, event() searches through the catalog of known sources
* to determine if a new source is more likely to appear first.
* Fourth:  if a new source is needed, we integrate over flux in m_currentRemainingFlux to find the 
* total, then integrate to a random number between 0 and that.  When the desired bin is reached, we subtract a single source from 
* m_currentRemainingFlux in the bin and declare a new source of that flux with addnewSource().
* Note:  it is possible, with this method, to have a fractional remaining source left in a bin.  When this source is called, the 
* program tries to compensate for it by removing a source from the adjacent lower bin as well.  This should be
* examined further.
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
#include <fstream>

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

      ///this should write a histogram of the constructed logN/logS characteristic 
      ///as defined by te existing sources.
      char* writeLogHistogram();

      /// write the characteristics of the current source distribution to a stream
      void writeSourceCharacteristic(std::ostream& out);

      ///remove a source from the list at the desired flux.
      void subtractFluxFromRemaining(double currentFlux);

      ///find the total remaining unresolved flux.
      void setFluxCharacteristics();
      
private:
    double m_totalFlux; // The total flux from the entire sky, particles/sec/m^2/steradian.
    double m_unclaimedFlux; // The amount of the flux "unaccounted for" by the known sources
    ///get a random flux from the remaining table of sources.
    double getRandomFlux();
    //the interpolated logN/logS characteristic
    long double logNlogS(long double flux);
    std::vector<std::pair<double,double> > m_logNLogS; // inputted logN.logS graph, to be used for finding random fluxes.
    std::vector<std::pair<double,double> > m_currentRemainingFlux;//remaining binned sources in the form (flux, number of sources)
    std::vector<PointSourceData> m_listOfDiffuseSources; //information on the sources which got added randomly.
    double m_minFlux; //the minimum flux to be considered in the spectrum
    double m_maxFlux; //the maximum flux to be considered.
    double m_totalIntegratedFlux; //integrated flux, using input file.
    std::string m_dtd; //the DTD document for the XML data to be read in.
    /// list of sources for easy lookup
    std::map<std::string, DOM_Element > m_sources;
};
