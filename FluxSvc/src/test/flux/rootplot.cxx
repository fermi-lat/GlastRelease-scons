// Flux test program that generates a ROOT macro to plot the flux
//

// $Header$

// Original author: Theodore Hierath


/**
  Test program for graphing the spectrums available through the flux
  package.
*/

#include "FluxSvc/FluxMgr.h"
#include "FluxSvc/EventSource.h"
#include "FluxSvc/ISpectrumFactory.h"
#include "FluxSvc/SpectrumFactoryTable.h"

#include "rootEnergyHist.h"
#include "rootAngleHist.h"

#include <fstream>

const int NUM_BINS = 30;
const int LOOP    = 30000;

const double ENERGY_MIN = 0.01;
const double ENERGY_MAX = 100.0;

static const char * default_arg="default";
static const char * default_graph="log";

void help() {
   std::cout << 
      "   Simple test program for a particle source.\n"
      "   Command line args are \n"
      "      '-sum' sums up the histograms'\n"
      "      '-bins <number of bins for histogram>'\n"
      "      '-events <number of events to create>'\n"
      "      '-min' or '-energy_min' sets minimum energy\n"
      "      '-max' or '-energy_max' sets maximum energy\n"
      "      '-list' lists the available spectra\n"
      "      '-file <file name>' writes energy vs. flux to filename.txt\n"
      "      '-trueflux' graphs using flux using per steradian\n"
      "      '-flux' graphs the flux vs. E instead of E*flux vs E\n"
      "      '-flux_min' manually set lower limit for graph\n"
      "      '-flux_max' manually set upper limit for graph\n"
      "      '-graph <log | semilogx | semilogy | linear>'\n"
      "      '-longsrc <sourcename>' for long-term energy averaging"
      "      '-help' for this help"
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
   DLL_DECL_SPECTRUM( FILESpectrum);
   
}

//void WARNING (const char * text ){  std::cerr << "WARNING: " << text << '\n';}
//void FATAL(const char* s){std::cerr << "\nERROR: "<< s;}

/**
  Test program for graphing the spectrums available through the flux
  package.

*/

int main(int argc, char** argv)
{
   int num_bins = NUM_BINS;  // Initialize to default number of bins
   int loop = LOOP;        // Initialize to default number of iterations
   int num_sources = 0;
   int num_longsources = 0;
   int current_arg = 1;
   int longtime=1;
   int longtimemax=20;
   double time=0.;  //time to use for flux and rate functions
   double energy_min = ENERGY_MIN;
   double energy_max = ENERGY_MAX;
   bool use_trueflux = false;
   bool use_flux = false;
   bool use_flux_min = false;
   bool use_flux_max = false;
   bool write_to_file = false;
   bool sum_mode = false;
   bool longterm=false;
   double flux_min;
   double flux_max;
   std::string arg_name(default_arg);
   std::string output_file_name;
   std::string ylabel_flux   = "Flux (particles/m^2/s/MeV/sr)";
   std::string ylabel_eflux  = "E*Flux (particles/m^2/s/sr)";
   std::string ylabel_sflux  = "Flux (particles/m^2/s/MeV)"; // Flux integrated over all solid angles
   std::string ylabel_seflux = "E*Flux (particles/m^2/s)";   // E*Flux integrated over all solid angles

   std::vector<std::string> sources;
   std::vector<int> longsources;
   
   flux_load();
   
   FluxMgr fm(sources); 

   
   // Process Command Line Arguments
   
   std::cout << "------------------------------------------------------" << std::endl;
   std::cout << " Flux test program: type 'rootplot -help' for help" << std::endl;
   std::cout << ( ( argc == 1)?  " No command line args, using defaults"
      :  "") << std::endl;
   
   
   
   while(current_arg < argc)
   {
      arg_name = argv[current_arg];
      if("-help" == arg_name || "help" == arg_name) 
      { 
         help(); 
         return 0; 
      }
      else if("-bins" == arg_name) 
         num_bins = atoi(argv[++current_arg]);
      else if("-events" == arg_name) 
         loop = atoi(argv[++current_arg]);
      else if("-min" == arg_name || "-energy_min" == arg_name) 
         energy_min = atof(argv[++current_arg]);
      else if("-max" == arg_name || "-energy_max" == arg_name) 
         energy_max = atof(argv[++current_arg]);
      else if("-flux_min" == arg_name)
      {
         use_flux_min = true;
         flux_min = atof(argv[++current_arg]);
      }
      else if("-flux_max" == arg_name)
      {
         use_flux_max = true;
         flux_max = atof(argv[++current_arg]);
      }
      else if("-flux" == arg_name)
         use_flux = true;
      else if("-trueflux" == arg_name)
         use_trueflux = true;
      else if("-file" == arg_name)
      {
         write_to_file = true;
         output_file_name = argv[++current_arg];
      }
      else if("-sum" == arg_name)
         sum_mode = true;
      else if("-list" == arg_name) 
      { 
         listSources(fm.sourceList());
         listSpectra(); 
         return 0; 
      }
      else if("-graph" == arg_name) 
         default_graph = argv[++current_arg];
      else if("-liny" == arg_name) default_graph = "semilogx";
      else if("-longsrc" == arg_name) {
         sources.push_back(argv[++current_arg]);
         //put the number of this source into the list for reference
         longsources.push_back(num_sources);
         num_sources++;
      }
      else if('-' == arg_name[0]) {std::cerr << "Unrecognized option "<< arg_name << ", -help for help" << std::endl;}
      else
      {
         sources.push_back(arg_name);
         num_sources++;
      }
      current_arg++;
   }
   
   
   // Use default source if no source was specified in arguments
   if(0 == sources.size())
   {
      sources.push_back(default_arg);
      num_sources++;
   }
   
   rootEnergyHist energy_hist(num_bins,energy_min,energy_max);
   rootAngleHist angle_hist(num_bins);
   
   // Process all the sources
   for(int i = 0; i < num_sources; i++)
   {
      int j;
      
      //decide whether or not this run should be longterm
      longterm = false;
      std::vector<int>::iterator longiter;
      for( longiter=longsources.begin(); longiter!=longsources.end() ;longiter++){
         if(*longiter==i) longterm=true;
      }

      if(longterm)
         fm.setExpansion(-1.);

          
     if((false == sum_mode && false==longterm)||(true==longterm && (longtime==1))) 
     {
        // Reset for new source
        energy_hist.reset();
        angle_hist.reset();
     }
     
     EventSource *e = fm.source(sources[i]);
     
     if(longterm){
        fm.pass(2.);
        time+=2.;
     }

     std::pair<double,double> loc=fm.location();
     std::cout << loc.first << "   " << loc.second << std::endl;
     //	  std::cout << "orbit angle=" << GPS::instance()-> << "orbit phase=" << << std::endl;
     
     if( 0==e ) {std::cerr << "Source \"" << sources[i] << "\" not found: -list for a list" << std::endl;
     return -1;}
     
     energy_hist.setGraphType(default_graph);
     energy_hist.setTitle( sources[i] );
     
     energy_hist.setXLabel("Kinetic Energy (GeV)");
     
     if(true == use_flux)
     {
        energy_hist.setFluxMode();
        if(false == use_trueflux)
           energy_hist.setYLabel(ylabel_sflux);
        else
           energy_hist.setYLabel(ylabel_flux);
     }
     else
     {
        if(false == use_trueflux)
           energy_hist.setYLabel(ylabel_seflux);
        else
           energy_hist.setYLabel(ylabel_eflux);
     }
     
     if(true == use_flux_min)
        energy_hist.setFluxMin(flux_min);
     
     if(true == use_flux_max)
        energy_hist.setFluxMax(flux_max);
     
     angle_hist.setGraphType(default_graph);
     angle_hist.setTitle( sources[i] );
     angle_hist.setPhiXLabel("Angle (degrees)");
     angle_hist.setPhiYLabel("Particles");
     angle_hist.setThetaXLabel("Cos(Theta)");
     angle_hist.setThetaYLabel("Particles");
     
     double scale_factor;

     if(false == use_trueflux)
     {
        double rate = e->rate(time)/e->totalArea();
        scale_factor = rate/loop*num_bins/log(energy_max/energy_min);
        std::cout << "RATE IS: " << rate << std::endl;
     }
     else
     {
        double flux = e->flux(time);
        scale_factor = flux/loop*num_bins/log(energy_max/energy_min);
        std::cout << "FLUX IS: " << flux << std::endl;
     }

     for(j = 0; j < loop; j++) 
     {
        FluxSource *f = e->event(time);
        double energy = f->energy();
        Vector dir = f->launchDir();
        double cos_theta = dir.z();
        
        double phi = atan2(dir.y(),dir.x());
        if(phi < 0)
           phi = 2*M_PI + phi;
        
        energy_hist.store(energy);
        angle_hist.storeTheta(cos_theta);
        angle_hist.storePhi(phi);
        
        if(j % 10000 == 0) std::cerr << "\r" << j << ": " << energy << "...";
     }

     std::cerr << "\n";
     // Scale factor must be applied before the draw member function can be called
     energy_hist.apply(scale_factor);
     angle_hist.apply(scale_factor);
     
     for(j = 0; j < num_bins; j++) 
     {
        std::cout << energy_min*pow(10,((j + 0.5)/num_bins) * energy_hist.retrieveRange() ) 
           << "   \t" << energy_hist.retrieveFlux(j) << "\t";
        std::cout << std::endl;         
     }
     
     
     if(false == sum_mode && false==longterm)
     {
        angle_hist.draw(1,"begin",i,num_sources);
        energy_hist.draw(1,"end",i,num_sources);

        delete e;
        std::cout << "Normal method" << std::endl;        
     }
     else if(true == sum_mode && i == num_sources - 1)
     {
        angle_hist.draw(1,"begin",0,1);
        energy_hist.draw(1,"end",0,1);

        delete e;
        std::cout << "Sum Mode method" << std::endl;
     }
     else if(longterm==true && longtime<=longtimemax)
     {
        longtime++;
        i--;
        std::cout << "Longterm run number " << longtime << std::endl; 
     }
     else
     {       
        if(true == write_to_file)
        {
           std::ofstream output_file;
           output_file_name += ".txt";
           output_file.open(output_file_name.c_str());
           
           energy_hist.writeFile(1./(longtime),output_file);
           output_file.close();
        }
        angle_hist.draw(1./double(longtime),"begin",i,num_sources);
        energy_hist.draw(1./double(longtime),"end",i,num_sources);

        delete e;
        longtime=1; 
               
     }

   } // for all sources   

   return 0;
}
