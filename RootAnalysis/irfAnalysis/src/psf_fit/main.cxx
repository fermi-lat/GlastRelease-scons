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

	//Perugia method to retreive the Xvalue of an histogram giving to it a y value

		Double_t GetXValue(TH1 *h,Double_t x_min,Double_t x_max)
{

  Int_t bin_min=h->GetXaxis()->FindBin(x_min);
  Int_t bin_max=h->GetXaxis()->FindBin(x_max);
  //cout<<"bins "<<bin_min<<' '<<bin_max<<endl; 
  Double_t Ymax=h->GetMaximum();
  Double_t diff=99999.,x_good=0.;
  for(int i=bin_min;i<=bin_max;i++){
   Double_t y_temp= h->GetBinContent(i);
   if(fabs(y_temp - (Ymax*0.33)) < diff){
     diff=fabs(y_temp - (Ymax*0.33));
     x_good= h->GetBinCenter(i);
   }
  }
   return x_good;

}

PsfModel(const std::string &outputFile, int maxTrys=5) 
        : Fitter(outputFile), m_maxTrys(maxTrys) 
    {

        m_func 
            = new TF1("myModel", 
                "x * ([0] * exp (- 0.5  * x * x /( [1] * [1]))  + [2] * exp (- (x * [3])^[4])   )"
            );
            
		m_gauss 
            = new TF1("myGaussModel", 
			"x * ([0] * exp (- 0.5  * x * x /( [1] * [1])) ) "
            );
			
		m_expo 
            = new TF1("myExpoModel", 
                "x * ( [2] * exp (- (x * [3])^[4])   )"
            );

		m_stgaus  
            = new TF1("myStdGauss", 
                "gaus"
            );
	
		//"x*( [0]*exp(-0.5*x*x/([1]*[1]))/[1] + [2]*exp(-x/[3]) )");  //+ [4]*exp(-x/[5]) )");

        // Set up the branches in the output tree.
        const char* names[]={"p0", "p1", "p2", "p3", "p4"};
        int npars=sizeof(names)/sizeof(void*);
        m_params.resize(npars);
        for( int i=0; i<npars; ++i){
            m_tree->Branch(names[i], &m_params[i], (std::string(names[i])+"/D").c_str());
        }
		

    }
    void applyFit(TH1 * h) 
    {
		// New fit from Perugia group with the pre-fit applied 

      Double_t histo_max=h->GetMaximum();
	  Double_t Xmax = h->GetXaxis()->GetXmax();
      Double_t Fitmax = Xmax*0.7; 
      Double_t Xgaus2 =0.;
      Double_t Xgaus1 =0.;
      Xgaus2= GetXValue(h,0.5,1.5);
      Double_t Xexpo1 =Xgaus2; 
      Double_t Xexpo2 =Xmax/2.;  
      Double_t par0gaus[3];
      Double_t par00gaus[3];
      Double_t parout[5];
 
//pre-fit of tails with an (x*exponential) function   
      Double_t par0expo[3]={histo_max*10.,40.,0.7}; 
      m_expo->SetParameters(par0expo[0],par0expo[1],par0expo[2]);
      m_expo->SetParLimits(0,histo_max,histo_max*100.);
      m_expo->SetParLimits(2,0.5,0.9);
      h->Fit("myExpoModel","","",Xexpo1,Xexpo2);
      m_expo->GetParameters(par0expo);

//pre-fit of peak with an (x*gaussian) function   
      h->Fit("myStdGauss","","",0.,Xgaus2);
      m_stgaus->GetParameters(par00gaus);
      m_gauss->SetParameters(histo_max,par00gaus[2]);
      m_gauss->SetParLimits(1,0.001,Xmax); //forced to be positive
      h->Fit("myGaussModel","","",Xgaus1,Xgaus2);
      m_gauss->GetParameters(par0gaus);
  
//Total fit
      m_func->SetParameters(par0gaus[0],par0gaus[1]*2.,par0expo[0],par0expo[1],par0expo[2]);
      m_func->SetParLimits(0,0.001,histo_max*10.);  //forced to be positive
      m_func->SetParLimits(2, 0, 10. );  //ditto
      m_func->SetParLimits(4,0.,1.5); 

	  h->Fit("myModel","","",0.,Fitmax);
    }
private:
   int m_maxTrys;
   TF1 * m_gauss;
   TF1 * m_expo;
   TF1 * m_stgaus;

};
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void fit_psfs();
void fit_edisp();

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
   
//   fit_edisp();   
//   fit_psfs();

   return 0;
}

void fit_edisp() {
   std::string gaussian("[0]/[2]*exp(-0.5*pow((x - [1])/[2], 2.))");
   double gauss_params[] = {1., 1., 0.2};
   
   std::vector<double> gaussParams(gauss_params, gauss_params+3);
   Fitter * myfit = new GenericFitter(gaussian, gaussParams, 
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
}

// void fit_psfs() {
//    std::string energyScaling("[0]*pow(1e-2/(1./McEnergy + [2]), [1])");

// //-------thin--------------
//    Fitter * myfit = new PsfModel("my_psf_thin_parameters.root");
//    MakeDists my_psf_thin("my_psf_thin.root");
//    my_psf_thin.set_user_cut(TCut("Tkr1FirstLayer<12"));

//    std::vector<double> params(3);
//    params[0] = 0.075;
//    params[1] = -0.85;
//    params[2] = 5e-4;
//    my_psf_thin.setEnergyScaling(energyScaling, params);

//    if (!my_psf_thin.fileExists()) {
//       my_psf_thin.project("BestDirErr", 0, 10, 100);
//    }
//    my_psf_thin.set_ymax(0.5); 
//    my_psf_thin.set_ymin(1e-4);
//    my_psf_thin.draw("my_psf_thin.ps", true , myfit);
//    delete myfit;

//    my_psf_thin.addCutInfo("my_psf_thin_parameters.root", "fitParams");
//    my_psf_thin.addEnergyScaling("my_psf_thin_parameters.root", 
//                                 "energyScaling");

// //-------thick-------------
//    myfit = new PsfModel("my_psf_thick_parameters.root");
//    MakeDists my_psf_thick("my_psf_thick.root");
//    my_psf_thick.set_user_cut(TCut("Tkr1FirstLayer>11"));
//    params[0] = 0.15;
//    params[1] = -0.85;
//    params[2] = 1e-4;
//    my_psf_thick.setEnergyScaling(energyScaling, params);

//    if (!my_psf_thick.fileExists()) {
//       my_psf_thick.project("BestDirErr", 0, 10, 100);
//    }
//    my_psf_thick.set_ymax(0.5); 
//    my_psf_thick.set_ymin(1e-4);
//    my_psf_thick.draw("my_psf_thick.ps", true , myfit);
//    delete myfit;

//    my_psf_thick.addCutInfo("my_psf_thick_parameters.root", "fitParams");
//    my_psf_thick.addEnergyScaling("my_psf_thick_parameters.root", 
//                                  "energyScaling");
// }
