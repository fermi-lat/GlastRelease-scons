// $Header$
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <utility>

// CLHEP
#include "CLHEP/Random/JamesRandom.h"

#include "../CrProton.h"
#include "../CrElectron.h"

const int BIN_NUM = 25;
const int LOOP    = 10000;

const double energy_min = 0.1;
const double energy_max = 100.0;

static const char * default_source="default";

void help() {
   std::cout << 
      "   Simple test program for a proton source.\n"
      "   Command line args are \n"
      "      'nbins <number of bins for histogram>'\n"
      "      'nevents <number of events to create>'\n"
      "      'graph' to generate a c++ script for ROOT\n"
      "      'help' for this help"
      << std::endl;
}

int main(int argc, char** argv)
{
   int bin_num = BIN_NUM;  // Initialize to default number of bins
   int loop = LOOP;        // Initialize to default number of iterations
   
   //CrProton src;
   CrElectron src;
   
   HepRandomEngine* engine = new HepJamesRandom;
   
   std::string source_name(default_source);
   
   if ( argc > 1 ) source_name = argv[1];
   if( source_name =="help") { help(); return 0; }
   if( source_name =="nbins") bin_num = atoi(argv[2]);
   if( source_name =="nevents") loop = atoi(argv[2]);
   if(argc > 3 ) source_name = argv[3];
   if( source_name =="nbins") bin_num = atoi(argv[4]);
   if( source_name =="nevents") loop = atoi(argv[4]);

   
   // Dynamically create and initialize bin array
   int *hist = new int[bin_num];
   for (int i = 0; i < bin_num; i++) hist[i] = 0;


   // Calculate constant logs
   double log10_energy_min = log10(energy_min);
   double log10_energy_max = log10(energy_max), range=log10_energy_max-log10_energy_min;

   double flux = src.flux(), scale_factor = flux/loop*bin_num/range;

   {for (int i = 0; i < loop; i++){
      src.selectComponent(engine);
      double energy = src.energySrc(engine); 
      std::pair<double,double> dir = src.dir(energy, engine);
      // Notice: dir consists of cos(zenith_angle) and azimuth [rad]
  
     hist[int(bin_num * (log10(energy) - log10_energy_min) / range)]++;


      if (i % 10000 == 0) std::cerr << '\015' << i << ": " << energy << "...";
   }}
   std::cerr << "\n";
   
   {for (int i = 0; i < bin_num; i++){
       std::cout << pow(10,((i + 0.5)/bin_num) * range + log10_energy_min) 
           << "   \t" << hist[i] << "\t";
       std::cout << std::endl;
   }}
   
   // Create file that can be run with ROOT to generate a graph

   std::ofstream out_file("graph.cxx");

   out_file << "{\n"
               "   gROOT->Reset();\n"
	       "   double emin= "<< energy_min << ", emax=" << energy_max << ";\n"; 
   out_file << "   double count[] = {\n";
   {for(int i = 0; i < bin_num; i++)
       out_file << std::setw(5) << hist[i] << (i%5==4? ",\n" : ",");
   }

   out_file << "      };\n\n"
               "   int bin_num = "<< bin_num << ";\n"
               "   double range = " << range << ";\n"
               "   double scale_factor =" << scale_factor << ", log10_energy_min=log10(emin);\n"
               "   double e_energy[bin_num];\n"
               "   double e_count[bin_num], energy[bin_num];\n\n"
               "   for(int i = 0; i < bin_num; i++) {\n"
               "      energy[i] = pow(10,((i+0.5)/bin_num) * range+ log10_energy_min);\n"
               "      e_count[i] = sqrt(count[i])*scale_factor;\n"
               "      e_energy[i] = 0;\n"
               "      count[i] *=scale_factor;\n"
               "   }\n"
               "   \n"
               "   //                 name, title, wtopx, wtopy, ww, wh\n"
               "   c1 = new TCanvas(\"c1\",\"Graph Window\", 200, 200, 700, 500);\n"
               "   c1->SetGrid();\n"
               "   c1->GetFrame()->SetFillColor(21);\n"
               "   c1->GetFrame()->SetBorderSize(12);\n"
               "   c1->SetLogx(1); c1->SetLogy(1);\n"
               "   gr = new TGraphErrors(bin_num,energy,count,e_energy,e_count);\n"
               "   gr->SetTitle(\"Proton Count vs. Energy\");\n"
               "   gr->SetMarkerColor(4);\n"
               "   gr->SetMarkerStyle(21);\n"
               "   gr->Draw(\"AP\");\n"
               "   TAxis* ax = gr->GetXaxis(), *ay= gr->GetYaxis();\n"
               "   ay->SetTitle(\"E*Flux (particles/m^2/s/)\"); ay->CenterTitle(1);\n"
	       "   ax->SetLimits(emin, emax);\n"
               "   ax->SetTitle(\"Kinetic Energy (GeV)\"); ax->CenterTitle(1); \n"
	       "   ax->SetTitleOffset(1.2);\n"
               "   c1->Modified();\n"
               "}\n";

   // Clean up
   
   out_file.close();
   system("root -l graph.cxx");

   delete engine;
   delete [] hist;
   
   return 0;
}
