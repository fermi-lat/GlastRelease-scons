#include <iterator>
#include <iostream>
#include <stdio.h>
#include <vector.h>

//include files for ROOT
#include "TROOT.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TH1.h"
#include "TRandom.h"

//GRB include files
#include "src/GRB/GRBShell.h"
#include "src/GRB/GRBShock.h"
#include "src/GRB/GRBConstants.h"
#include "src/GRB/GRBUniverse.h"
#include "src/GRB/GRBSim.h"

// ---------------- Main Pprogram ------

int main(int argc, char **argv){
  std::vector<double> spectrum;
  GRBConstants csta;
  GRBSim* _myGRB = new GRBSim;
  _myGRB->Start();
  
  // Initializing ROOT:
  //  TROOT root("GUI", "GUI test environement"); 
  //^ not need, inizialized in GRBSim!
  
  TApplication theApp("App", 0, 0);
  float phenergy;
  double enmin=1.0e+7;
  TH1D *hh1 = new TH1D("hh1","hh1",csta.enstep,log10(csta.enmin),log10(csta.enmax));
  TH1D *hh2 = new TH1D("hh2","hh2",csta.enstep,log10(csta.enmin),log10(csta.enmax));
  
  double ph;
  double eextra;
  double tim;
  int temp;
  double etot2;
  temp=1;
  
  TCanvas *cc1 = new TCanvas();  
  cc1->Divide(0,2);
  cout<<"Digit a time (in sec) or a 0 for the complete evolution:  "<<endl;
  while(temp==1){
    cin>>tim;
    cout<<"reset!!"<<endl;
    hh1->Reset();
    hh2->Reset();
    spectrum.clear();
    //for (int en=0;en<csta.enstep;en++){
    //      spectrum.push_back(0.0);
    //}
    
    if (tim>0){
      etot2=0.0;
      _myGRB->ComputeFlux(tim);
      spectrum=_myGRB->Spectrum();
      
      for (int en=0;en<csta.enstep;en++){
      	hh1->SetBinContent(en+1,spectrum[en]);
      }
      cout<<"Time  = "<< tim <<" Ftot  [eV/s/m^2]= "<<_myGRB->IFlux(csta.enmax/csta.enstep*en)<<endl;
      cout<<"                    Counts[ph/s/m^2]=" <<_myGRB->IRate(csta.enmax/csta.enstep*en)<<endl;
      cout<<"                    Energy[eV/m^2]=" <<_myGRB->IEnergy(csta.enmax/csta.enstep*en)<<endl;
      
      
      cc1->cd(1);
      hh1->Draw("AL");
      cc1->Update();      
      //      hh1->Integral(e1);
      eextra=0.0;
      int j=0;
      while (eextra<_myGRB->IEnergy(enmin)){
	phenergy=_myGRB->DrawPhotonFromSpectrum(spectrum,enmin);
	ph=log10(phenergy*1e+9);
	hh2->Fill(ph);
	eextra+=pow(10,ph);
	j++;
      }
   
    cout<< "Number of photons extractad: "<<j<<" Energy extracted = "<<eextra<<endl<<endl;
    cc1->cd(2);
    hh2->Draw();
    cc1->Update();
  } else {
    tim=-1.0;
      temp=0;
    }
  }
  cc1->Clear();
  cc1->Delete();
  cout<<temp<<endl;
  ///////////////////////////////////////////////////////
  //// if the argument is empty or different from zero, 
  /// the moovie is skipped !!
  _myGRB->TotalProperties(1);
  //////////////////////////////////////////////////////////
  _myGRB->Plot();
  theApp.Run();
  cin>>tim;
  delete _myGRB;
}


