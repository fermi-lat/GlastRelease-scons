// $Header$


#include "FluxSvc/FluxMgr.h"
#include "FluxSvc/FluxSource.h"
#include "FluxSvc/SpectrumFactoryTable.h"
#include "GPS.h"
#include "FluxException.h" // defines FATAL_MACRO

#include "dom/DOM_Document.hpp"
#include "dom/DOM_Element.hpp"
#include "xml/Dom.h"
#include "xml/IFile.h"

#include "Orbit.h"


#define DLL_DECL_SPECTRUM(x)   extern const ISpectrumFactory& x##Factory; x##Factory.addRef();



FluxMgr::FluxMgr(const std::vector<std::string>& fileList, std::string dtdname)
: m_dtd(dtdname.empty()? "$(FLUXSVCROOT)/xml/source.dtd" : dtdname )
{
    if( fileList.empty() ){
        defaultFile();
    }else{		
        init(fileList);		
    }	
}

void FluxMgr::defaultFile(){
    std::vector<std::string> input;
    // must find the source_library.xml file.
    // set up the xml document to use for initialization parameters
    static const char* initialization_document="source_library.xml";
    const char* flux_root = ::getenv("FLUXSVCROOT");
    std::string doc_path= (flux_root? std::string(flux_root)+"/xml/" : "");
    input.push_back(/*doc_path+initialization_document*/"$(FLUXSVCROOT)/xml/source_library.xml");	
    init(input);
}

void FluxMgr::init(const std::vector<std::string>& fileList){	
    std::string fileName;
    
    xml::XmlParser parser;
    std::vector<std::string>::const_iterator iter = fileList.begin();
    
    std::string xmlFileIn=  writeXmlFile( fileList);
    
    // a quick way of displaying what goes to the parser
    //std::cout << xmlFileIn <<std::endl;

    m_library_doc = parser.parse(xmlFileIn);
    
    if (m_library_doc == DOM_Document()) {
        FATAL_MACRO("Parse error: processing the document" << std::endl
            << xmlFileIn << std::endl);
        return;
    }
    
    // Root element is of type source_library.  Content is
    // one or more source elements.
    
    s_library = m_library_doc.getDocumentElement();
    
    // loop through the source elements to create a map of names, DOM_Elements
    if (s_library != DOM_Element()) {
        
        DOM_Element child = xml::Dom::getFirstChildElement(s_library);
        DOM_Element toplevel = xml::Dom::getFirstChildElement(s_library);
        while(child != DOM_Element()){
            char * b = xml::Dom::transToChar(child.getAttribute("name"));
            char * c = "";
            while(*b == *c){
                s_library = child;
                child = xml::Dom::getFirstChildElement(s_library);
                b = xml::Dom::transToChar(child.getAttribute("name"));
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
    // these are the spectra that we want to make available
    DLL_DECL_SPECTRUM( CHIMESpectrum);
    DLL_DECL_SPECTRUM( AlbedoPSpectrum);
    DLL_DECL_SPECTRUM( HeSpectrum);
    DLL_DECL_SPECTRUM( GalElSpectrum);
    DLL_DECL_SPECTRUM( CrElectron);
    DLL_DECL_SPECTRUM( CrProton);
    DLL_DECL_SPECTRUM( FILESpectrum);
    
}


FluxMgr::~FluxMgr(){
}

EventSource* FluxMgr::source(std::string name)
{
    // first check that it is in the library
    if( m_sources.find(name)==m_sources.end() ) {
        // nope. Maybe a Spectrum object
        Spectrum* s = SpectrumFactoryTable::instance()->instantiate(name);
        
        return s? new FluxSource(s) : (EventSource*)0;
    }
    return getSourceFromXML(m_sources[name]);
}


// sourceFromXML - create a new EventSource from a DOM element
// instantiated, e.g., from a description in source_library.xml
EventSource*  FluxMgr::getSourceFromXML(const DOM_Element& src)
{
    DOM_Node    childNode = src.getFirstChild();
    if (childNode == DOM_Node()) {
    /*
    FATAL_MACRO("Improperly formed XML event source");
    return 0;
        */
        // no child node: expect to find the name defined.
        return  new FluxSource(src);
    }
    
    DOM_Element sname = xml::Dom::getFirstChildElement(src);
    if (sname == DOM_Element() ) {
        FATAL_MACRO("Improperly formed XML event source");
        return 0;
    }
    // If we got here, should have legit child element
    if ((sname.getTagName()).equals("spectrum")) {
        
        return  new FluxSource(src);
    }
    else if ((sname.getTagName()).equals("nestedSource")) {
        
        // Search for and process immediate child elements.  All must
        // be of type "nestedSource".  There may be more than one.
        // Content model for nestedSource is EMPTY, so can omit check
        // for that in the code
        
        CompositeSource* cs = new CompositeSource();
        do { 
            //        DOM_Element sourceChild = (DOM_Element &) childNode;
            DOM_Element selem = 
                getLibrarySource(sname.getAttribute("sourceRef"));
            if (selem == DOM_Element()) {
                FATAL_MACRO("source name" << 
                    xml::Dom::transToChar(sname.getAttribute("sourceRef")) << 
                    "' not in source library");
            }
            cs->addSource(getSourceFromXML(selem));
            
            sname = xml::Dom::getSiblingElement(sname);
        } 
        while (sname != DOM_Element() );
        return cs;
    }
    else {
        FATAL_MACRO("Unexpected element: "<< 
            xml::Dom::transToChar(sname.getTagName()) );
    }
    return 0;
}



// source library lookup.  Each source is uniquely identified
// by its "name" attribute because "name" is of type ID
DOM_Element    FluxMgr::getLibrarySource(const DOMString& id)
{
    // quit if the library was unitialized
    if (s_library == DOM_Element() ) return DOM_Element(); 
    
    return m_library_doc.getElementById(id);
}

std::list<std::string> FluxMgr::sourceList() const
{
    std::list<std::string> s;
    for( std::map<std::string, DOM_Element>::const_iterator it = m_sources.begin();
    it != m_sources.end();
    ++it){
        s.push_back((*it).first);
    }
    return s;
}
/// generate some test output
void FluxMgr::test(std::ostream& cout, std::string source_name, int count)
{   
    EventSource* e = source(source_name);
    setExpansion(1.);
    double time=0.;

    cout << "running source: " << e->fullTitle() << std::endl;
    cout << " Total rate is: " << e->rate(0.) << " Hz into " << e->totalArea() << " m^2" << std::endl;
    //cout << "LaunchType" << f->retLaunch() << "Pointtype" << f->retPoint() <<std::endl;
    cout << "    Generating " << count << " trials " << std::endl;
    cout << " --------------------------------" << std::endl;
    
    //testing rotateangles function
    GPS::instance()->rotateAngles(std::make_pair<double,double>(0.0,0.3));

    FluxSource* f;
    for( int i = 0; i< count; ++i) {
        
        //testing - pass time
        pass(0.01);
        time+=0.01;

        f = e->event(time);
        //TESTING THE lat, lon FUNCTIONS
        //cout << std::endl << "lat=" << GPS::instance()->lat() << ' ' <<"lon=" << GPS::instance()->lon() << std::endl;
        //double curTime=GPS::instance()->time();
        //cout << std::endl << "testlat=" << GPS::instance()->orbit()->testLatitude(curTime) << ' ' << "testlon=" << GPS::instance()->orbit()->testLongitude(curTime) << std::endl;

        cout << "LaunchType = " << f->refLaunch() << " , Pointtype = " << f->refPoint() <<std::endl;
        cout << f->spectrum()->particleName();
        cout << "(" << f->energy();
        cout << " GeV), Launch: " << f->launchPoint() 
            << " Dir " << f->launchDir() << " ,Flux="
            << f->flux(time) << " ,Interval="
            << f->interval(time) << std::endl;
    }
    cout << "------------------------------------------------------" << std::endl;
    
    delete e;
    
}


void FluxMgr::addFactory(std::string name, const ISpectrumFactory* factory ) {
    SpectrumFactoryTable::instance()->addFactory(name,factory);
}

void FluxMgr::setGlastAngles(std::pair<double,double> ang){
    GPS::instance()->rotateAngles(ang);
}


void FluxMgr::setGlastPosition(std::pair<double,double> pos){
    GPS::instance()->ascendingLon(pos.first);
    GPS::instance()->ascendingLon(pos.second);
}

void FluxMgr::setExpansion (double p){
    // set the expansion factor for the orbit (-1) = random
    GPS::instance()->expansion(p);
}

// pass a specific amount of time
void FluxMgr::pass(double t){
    GPS::instance()->pass(t);
    synch();
}

GPStime FluxMgr::time () const{
    return GPS::instance()->time();
}

void FluxMgr::synch(){
    GPS::instance()->synch();
}

void sampleintvl ( /*GPStime*/double t ){
    GPS::instance()->sampleintvl(t);
}

//get the current satellite location
std::pair<double,double> FluxMgr::location(){
    return std::make_pair<double,double>(GPS::instance()->lat(),GPS::instance()->lon());
}

//get the transformation matrix.
Rotation FluxMgr::CELTransform(double time){
    return GPS::instance()->orbit()->CELtransform(time);
}


/** creates a document of the form

<?xml version='1.0' ?>
<!DOCTYPE source_library SYSTEM "d:\users\burnett\pdr_v7r1c\flux\v5r3/xml/source.dtd" [
<!ENTITY librarya SYSTEM "d:\users\burnett\pdr_v7r1c\flux\v5r3/xml/user_library.xml" >
<!ENTITY libraryb SYSTEM "d:\users\burnett\pdr_v7r1c\flux\v5r3/xml/source_library.xml" >
]>
<source_library>
&librarya;
&libraryb;
</source_library>

*/
std::string FluxMgr::writeXmlFile(const std::vector<std::string>& fileList) {
  
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
    fileString << "<?xml version='1.0' ?>" << std::endl << "<!DOCTYPE source_library" 
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
    
    fileString << "]>" << std::endl << "<source_library>" << std::endl;
    iter = fileList.begin();
    libchar = 'a';

    //as long as there are files in the file list...
    for (;iter != fileList.end(); iter++) {
        // add a reference to the file name
        fileString << "&library" << libchar << ";" << std::endl;       
        libchar++;
    }
    
    fileString << "</source_library>" << '\0';
    return fileString.str();

}
