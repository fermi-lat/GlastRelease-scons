#include "CompositeDiffuse.h"
#include "CLHEP/Random/RandFlat.h"
#include "FluxSvc/FluxSource.h"
#include "SimpleSpectrum.h"
#include <strstream>

void CompositeDiffuse::addSource (EventSource* aSource)
{
    m_sourceList.push_back(aSource);
    //flux should ensure proper units.
    // want EVERYTHING in particles/sec/m^2/steradian.
    double flux = aSource->flux(EventSource::time());
    //std::cout << "adding source of flux " << flux << std::endl;
    //EventSource::setFlux( flux );
    m_unclaimedFlux-=flux;
}

FluxSource* CompositeDiffuse::event (double time)
{
    //Purpose: To determine which contained source is currently represented
    //Input: Current time
    //Output: Pointer to the current represented FluxSource object.
    EventSource::setTime(time);
    
    int m_numofiters=0;
    
    // here we should be setting the total rate as the maximum sum rate of everything - FIX!
    // double mr = rate(EventSource::time());
    double mr = m_totalFlux;
    
    // do this once if there is no source, or no rate at all (null source?)
    if( m_sourceList.size()==0 || mr == 0) {
        m_recent = m_sourceList.front();
    }else {
        
        // more than one:: choose on basis of relative rates
        double  x = RandFlat::shoot(mr), y = 0;
        std::vector<EventSource*>::iterator  now = m_sourceList.begin();
        std::vector<EventSource*>::iterator  it = now;
        
        double intrval=0.,intrmin=100000.;
        
        for ( int q=0; now != m_sourceList.end(); ++now) {
            (*now)->event(time); // to initialize particles, so that the real interval for the particle is gotten.
            intrval=(*now)->interval(EventSource::time());
            
            if(intrval < intrmin){
                it=now;
                intrmin=intrval;
                m_numofiters=q;
            }
            
            m_recent = (*it);
            q++;
        }
        //so now m_recent is the "soonest" source, and intrmin is the "soonest": time.
        // but, what if the "leftover" flux sends a source even faster?
        intrval=remainingFluxInterval();
        if(intrval < intrmin){
            
            addNewSource();
            intrmin=intrval;
            //++now;
            //m_recent = *(m_sourceList.end());
            it = m_sourceList.end();
            --it;
            m_recent = (*it);
        }
        //set the interval that won out after everything.
        setInterval(intrmin);
    }
    //update the time
    //m_time += interval(m_time);
    // now ask the chosen one to generate the event, if there is a rate
    return (FluxSource*)m_recent;//->event(time);
}


void CompositeDiffuse::addNewSource(){
    //here, set the flux in a random fashion, determined by the log N/ log S characteristic graph
    double flux=1.0;
    flux = getRandomFlux();
    //Set up the spectrum with the appropriate power spectrum
    Spectrum * spec=new SimpleSpectrum("gamma", 100.0);
    
    // Now make the random vector correesponding to random galactic direction
    double  costh = -RandFlat::shoot(-1.,1./*_minCos, _maxCos*/);
    double  sinth = sqrt(1.-costh*costh);
    double  l = RandFlat::shoot(-180, 180);       
    double b = (acos(costh)*180/M_PI)-90.;
    EventSource* aSource=new FluxSource(flux,spec,l,b);
    PointSourceData thissource; 
    thissource.flux = flux;
    thissource.l = l;
    thissource.b = b;
    //thissource.energyIndex = spec;
    m_listOfDiffuseSources.push_back(thissource);
    
    FluxSource* abc = (FluxSource*)aSource;
    //then add it into the list of sources..
    addSource(aSource);
    //..and subtract the total flux from what remains...
    m_unclaimedFlux -= flux;
}


double CompositeDiffuse::remainingFluxInterval(){
    //std::cout << "unclaimed flux is " << m_unclaimedFlux << std::endl;
    
    double  r = (solidAngle()*(m_unclaimedFlux)*6.);
    
    if (r <= 0){ return 1E20;
    }else{  
        double p = RandFlat::shoot(1.);
        return (-1.)*(log(1.-p))/r;
    }
    
}

long double CompositeDiffuse::logNlogS(long double flux){
    //READ MY LIPS - THIS NEEDS TO DO A PROPER INTERPOLATION!!!
    //Purpose:  this is designed to interpolate over the logN/logS characteristic.
    //Input:  the flux (in linear units)
    //Output:  the number of sources per 1/5 decade, at the input flux.
    //Caveats: the vector should be sorted - right now the file needs to be in order of ascending flux
    int i = 0;
    double logHighFlux,logLowFlux;
    i++; //there needs to be at least a space of 1 data point for measurement.
    while( (log10(flux) - m_logNLogS[i].first >= 0) && (m_logNLogS[i].first!=m_maxFlux) ){
        i++;
    }
    logHighFlux = m_logNLogS[i].first;
    double logHighSources = m_logNLogS[i].second;
    logLowFlux = m_logNLogS[i-1].first;
    double logLowSources = m_logNLogS[i-1].second;
    
    double logNumSources = logLowSources + ((logHighSources-logLowSources)*(log10(flux)-logLowFlux));
    //std::cout << "loglowsources = " << logLowSources << std::endl;
    
    return pow(10,logNumSources);
}

double sizeOfFifthDecade(double currentFlux){
    return /*currentFlux**/(pow(10,0.2)-1.);
}

double CompositeDiffuse::getRandomFlux(){
    //NEED TO SET TOTAL INTEGRATED FLUX!!!!
    long double prob=RandFlat::shoot(m_totalIntegratedFlux);
    long double dx=0.0000000000001;
    
    double currentFlux = m_minFlux;
    while(prob > 0 && currentFlux<m_maxFlux){
        //dx=0.01*pow(10.0,log10(i));
        double sourcesPerFlux = logNlogS(currentFlux)/sizeOfFifthDecade(currentFlux);
        prob-=(dx)* sourcesPerFlux;
        currentFlux+=sourcesPerFlux*dx;
        //	printf("\nin the findandaddnew loop; i=%12.10e, prob=%12.10e, dx=%12.10e , logi=%lf\n",i,prob,dx,log10(i));
    }
    std::cout << "New Source created in CompositeDiffuse, with flux = " << currentFlux << std::endl;
    return currentFlux;
}

void CompositeDiffuse::setFileFlux(){
    double  currentFlux = m_minFlux ;
    m_totalIntegratedFlux = 0.;  //to initialize
    long double dx=0.0000000000001; 
    while(currentFlux < m_maxFlux){
        double sourcesPerFlux = logNlogS(currentFlux)/sizeOfFifthDecade(currentFlux);
        m_totalIntegratedFlux+=(dx)* sourcesPerFlux;
        currentFlux+=sourcesPerFlux*dx;
        //std::cout << "m_currentFlux = " << currentFlux << std::endl;
        
    }
    
}

std::string CompositeDiffuse::writeXmlFile(const std::vector<std::string>& fileList) {
/** purpose: creates a document of the form

  <?xml version='1.0' ?>
  <!DOCTYPE data_library SYSTEM "d:\users\burnett\pdr_v7r1c\flux\v5r3/xml/source.dtd" [
  <!ENTITY librarya SYSTEM "d:\users\burnett\pdr_v7r1c\flux\v5r3/xml/user_library.xml" >
  <!ENTITY libraryb SYSTEM "d:\users\burnett\pdr_v7r1c\flux\v5r3/xml/source_library.xml" >
  ]>
  <data_library>
  &librarya;
  &libraryb;
  </data_library>
  
    */
    std::strstream fileString;
    // Unique tag to add to ENTITY elements in the DTD.
    char libchar = 'a';
    std::string inFileName;
    
    std::vector<std::string>::const_iterator iter = fileList.begin();
    
    //the default DTD file
    inFileName=m_dtd;
    //replace $(FLUXROOT) by its system variable
    xml::IFile::extractEnvVar(&inFileName);
    
    //this stuff goes in the beginnning of the XML file to be read into the parser
    fileString << "<?xml version='1.0' ?>" << std::endl << "<!DOCTYPE data_library" 
        << " SYSTEM " << '"' << inFileName << '"' << " [" << std::endl;
    
    //as long as there are files in the file list...
    for (;iter != fileList.end(); iter++) {
        
        // get the file name, and evaluate any system variables in it
        inFileName=(*iter).c_str();
        xml::IFile::extractEnvVar(&inFileName);
        
        //then make an ENTITY entry specifying where the file is
        fileString << "<!ENTITY " << "library" << libchar << " SYSTEM " << '"' 
            << inFileName << "\" >" << std::endl;      
        libchar++;
    }
    
    fileString << "]>" << std::endl << "<data_library>" << std::endl;
    iter = fileList.begin();
    libchar = 'a';
    
    //as long as there are files in the file list...
    for (;iter != fileList.end(); iter++) {
        // add a reference to the file name
        fileString << "&library" << libchar << ";" << std::endl;       
        libchar++;
    }
    
    fileString << "</data_library>" << '\0';
    return fileString.str();
    
}

void CompositeDiffuse::fillTable(){
    
    m_dtd = "$(FLUXSVCROOT)/xml/FluxSvcParams.dtd";
    std::string xmlFile = "$(FLUXSVCROOT)/xml/FluxSvcParams.xml";
    
    xml::XmlParser parser;
    
    //prepare the XML file
    std::vector<std::string> fileList;
    fileList.push_back(xmlFile);
    
    //then create the parser input.
    std::string xmlFileIn = writeXmlFile( fileList);
    
    // a quick way of displaying what goes to the parser
    //std::cout << xmlFileIn <<std::endl;
    
    DOM_Document m_library_doc = parser.parse(xmlFileIn);
    
    if (m_library_doc == DOM_Document()) {
        std::cout << "Parse error: processing the document" << std::endl
            << xmlFileIn << std::endl;
        return;
    }
    
    // Root element is of type data_library.  Content is
    // one or more datatable elements.
    
    DOM_Element s_library = m_library_doc.getDocumentElement();
    
    // loop through the data elements to create a map of names, DOM_Elements
    if (s_library != DOM_Element()) {
        
        DOM_Element child = xml::Dom::getFirstChildElement(s_library);
        DOM_Element toplevel = xml::Dom::getFirstChildElement(s_library);
        
        while (child != DOM_Element()) {
            while (child.getAttribute("name") == DOMString()) {
                s_library = child;
                child = xml::Dom::getFirstChildElement(s_library);
            }
            
            while (child != DOM_Element()) {
                std::string name = xml::Dom::transToChar(child.getAttribute("name"));
                m_sources[name]=child;
                child = xml::Dom::getSiblingElement(child);
            }
            
            child = xml::Dom::getSiblingElement(toplevel);
            toplevel=child;
        }
        
    }
    
    //now get the desired data table:
    DOM_Element characteristic = m_sources["logNlogScharacteristic"];
    
    
    DOM_Element sname = xml::Dom::getFirstChildElement(characteristic);
    if (sname == DOM_Element() ) {
        std::cout << "Improperly formed XML event source" << std::endl;
        return;
    }
    
    DOM_Element toplevel=sname;
    // If we got here, should have legit child element
    while (/*(sname.getTagName()).equals("logNlogSdata")*/sname != DOM_Element()) {
        double logflux = atof(xml::Dom::transToChar(sname.getAttribute("logflux")));
        double lognumsources = atof(xml::Dom::transToChar(sname.getAttribute("lognumsources")));
        m_logNLogS.push_back(std::make_pair<double,double>(logflux, lognumsources));
        
        //std::cout << logflux << " , " << lognumsources << std::endl;
        sname = xml::Dom::getSiblingElement(toplevel);
        toplevel=sname;
    }
    //now figure out the maximum and minimum fluxes.
    m_minFlux = pow(10,(*m_logNLogS.begin()).first);
    m_maxFlux = pow(10,(m_logNLogS.back()).first);
    //std::cout << m_minFlux << "   " << m_maxFlux<<std::endl;
    return;          
}


///this should write a histogram of the constructed logN/logS characteristic 
///as defined by te existing sources.
char* CompositeDiffuse::writeLogHistogram(){
    
    std::vector<std::pair<double,double> > currentHistPoints;
    //the histogram starts at minFlux, ends at maxFlux, and holds number of events.
    for(double curFlux = m_minFlux ; curFlux*(1.+sizeOfFifthDecade(curFlux)) <=m_maxFlux ; curFlux+=curFlux*sizeOfFifthDecade(curFlux) ){
        currentHistPoints.push_back(std::make_pair<double,double>(curFlux,0.));  
    }
    
    std::vector<PointSourceData>::iterator now= m_listOfDiffuseSources.begin();
    for ( ; now != m_listOfDiffuseSources.end(); now++) {
        //std::vector<std::pair<double,double> >::iterator iter;
        for(int i=0 ; i <= currentHistPoints.size() ; i++){
            if(currentHistPoints[i].first<=(*now).flux && currentHistPoints[i+1].first>=(*now).flux) currentHistPoints[i].second++;
        }
    }
    //then we want to change the histogram bins to be centered in their ranges, and log everything
    std::vector<std::pair<double,double> >::iterator iter;
    for(iter=currentHistPoints.begin() ; iter!=currentHistPoints.end() ; iter++){
        (*iter).first = log10((*iter).first);
        (*iter).second = log10((*iter).second);
        (*iter).first += 1/10.;
    }
    //now, the output:
    std::strstream out2;// = new std::strstream;
    for(iter=currentHistPoints.begin() ; iter!=currentHistPoints.end() ; iter++){
        out2 << (*iter).first << '\t' << (*iter).second << std::endl;
    }
    std::cout << out2.str(); //WHY IS THERE JUNK ON THE END OF IT?
    return out2.str();
    
}


/// write the characteristics of the current source distribution to a stream
void CompositeDiffuse::writeSourceCharacteristic(std::ostream& out){
    out << fullTitle() << std::endl;
    out << "Diffuse Added Point Sources:" << std::endl;
    out << "l" << '\t' << "b" <<'\t' << "flux" << '\t'<< "energyIndex" << std::endl << std::endl;
    
    std::vector<PointSourceData>::iterator now= m_listOfDiffuseSources.begin();
    for ( ; now != m_listOfDiffuseSources.end(); now++) {
        out << (*now).l << '\t' <<(*now).b << '\t' << (*now).flux << '\t' << (*now).energyIndex << std::endl;
    }
    out << "histogram of reconstructed logN/logS:" << std::endl;
    out << "flux(integrated over 1/5 decade)" << '\t' << "number of sources" << std::endl << std::endl;
    out << writeLogHistogram() << std::endl;
    
}