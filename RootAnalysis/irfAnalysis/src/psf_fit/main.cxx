/** @file main.h
$Header$
*/

#include "MakeDists.h"
#include "Fitter.h"
#include "GenericFitter.h"
#include "TF1.h"
#include <cmath> // for mod

class PsfModel : public Fitter {
public:

   PsfModel(const std::string &outputFile, int maxTrys=5) 
      : Fitter(outputFile), m_maxTrys(maxTrys) {
      std::string psfFunc
         = "x*( [0]*exp(-0.5*x*x/[1]/[1]) + [2]*exp(-pow(x*[3], [4])) )";
      m_func = new TF1("myModel", psfFunc.c_str());

// Set up the branches in the output tree.
      const char* names[]={"p0", "p1", "p2", "p3", "p4"};;
      int npars=sizeof(names)/sizeof(void*);
      m_params.resize(npars);
      for (int i=0; i<npars; ++i) {
         m_tree->Branch(names[i], &m_params[i], 
                        (std::string(names[i])+"/D").c_str());
      }
      m_func->SetParLimits(0, 0.1, 1);  // minimum for the gaussian component
      m_func->SetParLimits(1, 0.2, 5);
      m_func->SetParLimits(2, 0, 10);;
      m_func->SetParLimits(3, 0, 10);;
      m_func->SetParLimits(4, 0.5, 1.5);;

   }
   ~PsfModel(){delete m_func;}

   void applyFit(TH1 * h) {
      // These initial parameters were obtained interactively
      double pinit[]={0.01, 1.0, 1.0, 1.0, 1.0};
      std::copy(pinit, pinit+sizeof(pinit)/sizeof(double), m_params.begin());
      m_func->SetParameters(&m_params[0]);
      int fitTrys = 0;
      while (h->Fit("myModel") && fitTrys++ < m_maxTrys);
   }

private:
    int m_maxTrys;
};
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
int main(){

   //-------thin --------------
   Fitter * myfit = new PsfModel("psf_thin_parameters.root");
   MakeDists psf_thin("psf_thin.root");
   psf_thin.set_user_cut(TCut("Tkr1FirstLayer<12"));
   if ( !psf_thin.fileExists() )
      psf_thin.project("BestDirErr/PSFscaleFactor", 0, 10, 100 );
   psf_thin.set_ymax(0.2); psf_thin.set_ymin(1e-4);
   psf_thin.draw("psf_thin.ps", true , myfit);
   delete myfit;
   psf_thin.addCutInfo("psf_thin_parameters.root", "fitParams");
   std::vector<double> Perugia_thin(1);
   Perugia_thin[0] = 0.08;
   psf_thin.setEnergyScaling("", Perugia_thin, true);
   psf_thin.addEnergyScaling("psf_thin_parameters.root", "energyScaling");
                             
   //-------thick --------------
   myfit = new PsfModel("psf_thick_parameters.root");
   MakeDists psf_thick("psf_thick.root");
   psf_thick.set_user_cut(TCut("Tkr1FirstLayer>11"));
   if ( !psf_thick.fileExists() )
      psf_thick.project("BestDirErr/PSFscaleFactor", 0, 10, 100 );
   psf_thick.set_ymax(0.5); psf_thick.set_ymin(1e-4);
   psf_thick.draw("psf_thick.ps", true , myfit);
   delete myfit;
   psf_thick.addCutInfo("psf_thick_parameters.root", "fitParams");
   std::vector<double> Perugia_thick(1);
   Perugia_thick[0] = 0.14;
   psf_thick.setEnergyScaling("", Perugia_thick, false);
   psf_thick.addEnergyScaling("psf_thick_parameters.root", "energyScaling");
   
   //-------energy dispersion-------
   std::string gaussian("[0]/[2]*exp(-0.5*pow((x - [1])/[2], 2.))");
   double gauss_params[] = {1., 1., 0.2};
   
   std::vector<double> gaussParams(gauss_params, gauss_params+3);
   myfit = new GenericFitter(gaussian, gaussParams, 
                             "edisp_thin_parameters.root");
   MakeDists edisp_thin("edisp_thin.root");
   edisp_thin.set_user_cut(TCut("Tkr1FirstLayer<12"));
   if (!edisp_thin.fileExists()) 
      edisp_thin.project("EvtEnergySumOpt/McEnergy", -0.5, 2.5, 100);
   edisp_thin.set_ymax(0.5);
   edisp_thin.set_ymin(1e-4);
   edisp_thin.draw("edisp_thin.ps", true, myfit);
   delete myfit;
   edisp_thin.addCutInfo("edisp_thin_parameters.root", "fitParams");
   
   gaussParams[0] = 1.;
   gaussParams[1] = 1.;
   gaussParams[2] = 0.2;
   myfit = new GenericFitter(gaussian, gaussParams, 
                             "edisp_thick_parameters.root");
   MakeDists edisp_thick("edisp_thick.root");
   edisp_thick.set_user_cut(TCut("Tkr1FirstLayer>11"));
   if (!edisp_thick.fileExists()) 
      edisp_thick.project("EvtEnergySumOpt/McEnergy", -0.5, 2.5, 100);
   edisp_thick.set_ymax(0.5);
   edisp_thick.set_ymin(1e-4);
   edisp_thick.draw("edisp_thick.ps", true, myfit);
   delete myfit;
   edisp_thick.addCutInfo("edisp_thick_parameters.root", "fitParams");
   
   return 0;
}

