// $Header$

#include "../../FluxMgr.h"
#include "FluxSvc/EventSource.h"
#include "FluxSvc/ISpectrumFactory.h"
#include "FluxSvc/SpectrumFactoryTable.h"
#include <iostream>
#include <algorithm>
#include "../../CHIMESpectrum.h"
#include "../../Orbit.h"


static int default_count = 200 ;
//Testing
static const char * default_source="default";
//Default
//static const char * default_source="CrElectron";

void help() {
    std::cout << 
        "   Simple test program to create a particle source, then run it.\n"
        "   Command line args are \n"
        "      <source name> <count>\n"
        "   with defaults \"" 
        <<  default_source << "\"," << default_count
        << "\n  Also, 'help' for this help, 'list' for a list of sources and Spectrum objects "
        << std::endl;
}
void listSources(const std::list<std::string>& source_list ) {
    std::cout << "List of available sources:" << std::endl;
    for( std::list<std::string>::const_iterator it = source_list.begin(); 
    it != source_list.end(); ++it) { 
        std::cout << '\t'<< *it << std::endl;
    }
    
}
void listSpectra() {
    std::cout << "List of loaded Spectrum objects: " << std::endl;
    std::list<std::string> spectra(SpectrumFactoryTable::instance()->spectrumList());
    for( std::list<std::string>::const_iterator it = spectra.begin(); 
    it != spectra.end(); ++it) { 
        std::cout << '\t'<< *it << std::endl;
    }
}


#define DLL_DECL_SPECTRUM(x)   extern const ISpectrumFactory& x##Factory; x##Factory.addRef();

void flux_load() {
  
  // these are the spectra that we want to make available
  DLL_DECL_SPECTRUM( CHIMESpectrum);
  DLL_DECL_SPECTRUM( AlbedoPSpectrum);
  DLL_DECL_SPECTRUM( HeSpectrum);
  DLL_DECL_SPECTRUM( GalElSpectrum);
  DLL_DECL_SPECTRUM( CrElectron);
  DLL_DECL_SPECTRUM( CrProton);
  //  DLL_DECL_SPECTRUM( GRBSpectrum);
  //  DLL_DECL_SPECTRUM( CREMESpectrum);
}

int main(int argn, char * argc[]) {
  using std::cout;
    using std::endl;
    flux_load();
    //testing 
    std::string test("23class");
    int ss = test.find_first_not_of("0123456789");
    
    int count = default_count;
    std::string source_name(default_source);
    
    
    //TESTING MULTIPLE XML INPUT
    std::vector<std::string> fileList;
    fileList.push_back("$(FLUXSVCROOT)/xml/user_library.xml");
    fileList.push_back("$(FLUXSVCROOT)/xml/source_library.xml");	
    FluxMgr fm(fileList);
    
    //FluxMgr fm;
    
    //Testing the addfactory function
    //    static PencilBeam* sean=TestSpec::instance();
    //    fm.addFactory("seantest", sean );
    //End Test
    
    if ( argn >1 ) source_name = argc[1];
    if( source_name =="help") { help(); return 0; }
    if( source_name =="list") { 
        listSources(fm.sourceList());
        listSpectra(); return 0; }
    if ( argn >2 ) count = ::atoi(argc[2]);
    
    cout << "------------------------------------------------------" << endl;
    cout << " Flux test program: type 'help' for help" << endl;
    cout << ( ( argn ==1)?  " No command line args, using default flux \""
        :  " Selected source name \"");
    cout  << source_name <<"\"" << endl;
    
    
    std::list<std::string> source_list(fm.sourceList());
    
    if( std::find(source_list.begin(), source_list.end(), source_name)==source_list.end() ) {
        std::list<std::string> spectra(SpectrumFactoryTable::instance()->spectrumList());
        
        if( std::find(spectra.begin(), spectra.end(), source_name)==spectra.end() ) {
            std::cout << "Source \"" << source_name << "\" not found in the list or sources!" << std::endl;
            listSources(source_list);
            std::cout << "or in specra list, which is:\n";
            listSpectra();
            
            return -1;
        }
    }
    // now have FluxMgr create and run it.
    fm.test(std::cout, source_name, count);
    return 0;
    
}

