#include <iterator>
#include <iostream>
#include <stdio.h>
#include <algorithm>
#include <vector>
// include files for ROOT
#include "TROOT.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
// my include files
#include "GRBShell.h"
#include "GRBShock.h"
#include "GRBConstants.h"
#include "GRBUniverse.h"
#include "GRBSim.h"


/*------------------------------------------------------*/
class ShockCmp{
public:
  bool operator()(GRBShock* Sho1, GRBShock* Sho2)
  {
    return Sho1->tobs()<Sho2->tobs();    
  }
};
/*------------------------------------------------------*/


GRBSim::GRBSim(){
  cout<<"******Staring The GRB Simulation******"<<endl;
  cout<<"******     Initializing ROOT    ******"<<endl;
  static const   TROOT root("GRB", "GRB Simulator");
  myUni=new GRBUniverse();
  cout<<"Dist  of the source = "<<myUni->Dist();
  double temp1=(cst.enmax/cst.enmin);
  double temp2=(1.0/cst.enstep);
  double denergy = pow(temp1,temp2);
  
  for(int en=0;en<=cst.enstep;en++)
    {
      m_energy[en]=cst.enmin*pow(denergy,en); 
      m_gammae[en]=(1.+cst.enmin*pow(denergy,en)/(1e+6*cst.mec2));
    }
  for(int en=0;en<cst.enstep;en++)
    {
      m_de[en]= m_energy[en+1]-m_energy[en];
    }
}
/*------------------------------------------------------*/
GRBSim::~GRBSim(){
  cout<<"*******Exiting The GRB Simulation"<<endl;
}
/*------------------------------------------------------*/
void GRBSim::Start() {
   /// Step 1: Creation of the shells:
  double ei = cst.etot/cst.nshell; //erg
  double dr=(cst.r0*cst.dr0);
  for(int i=cst.nshell;i>0;i--) {
    GRBShell* iShell = new GRBShell(ei);

    iShell->setThickness(dr);
    iShell->setRadius(i*cst.r0+dr);

    cout << " Shell n: "<<cst.nshell-(i-1)<<" Gamma= "<<iShell->Gamma()<<" Initiaml Radius "<<iShell->Radius()<< " Initial Thickness= " <<iShell->Thickness()<< endl;
    theShells.push_back(iShell);
  }
  double ssum;
  double tmax = cst.dt1*cst.nstep; 
  double time = 0.0;
  int nshock=0;

  /// Step 2: Calculation of the evolution:
  while(time<tmax)
    {
      for(int i=1;i<cst.nshell-nshock;i++){
	theShells[i]->evolve(time);      
      }
      for(int i=2;i<cst.nshell-nshock;i++) 
	{
	  GRBShell* Sh1=theShells[i-1];
	  GRBShell* Sh2=theShells[i];
	  if(Sh1->Radius()<=(Sh2->Radius()+Sh2->Thickness())) 
	    {
	      GRBShock* iShock = new GRBShock(Sh1, Sh2);
	      iShock->setTime(time);
	      theShocks.push_back(iShock);	
	      theShells.erase(&theShells[i]);
	      nshock++;
	    }
	}
      time+=cst.dt1;
    }
  cout<< "Number of Shocks = " <<nshock<< endl;
  if (nshock==0){
    cout<< "Sorry no shocks events!! "<< endl;
    std::exit(1);
  }
  
/*------------------------------------------------------*/
  /// Step 3: Sorting the shocks and setting t min=0
  std::sort(theShocks.begin(), theShocks.end(), ShockCmp());
  double t0 = theShocks[0]->tobs();
  for(int i=0;i<nshock;i++)
    {
      double temp=theShocks[i]->tobs();
      theShocks[i]->setTobs(temp-t0);
    cout<<"Shock num " << i << " @ time obs = " << theShocks[i]->tobs()<<endl;
    }
  // attention tmax redeclared here!!
  m_tmax=1.2*theShocks[nshock-1]->tobs()+0.1;
  m_dt=m_tmax/cst.nstep;
  cout<<cst.flagIC<<endl;
  for(int i=0;i<nshock;i++)
    {
      ssum=0.0;
      for (int tt=0;tt<cst.nstep;tt++){
	for (int en=0;en<cst.enstep;en++){
	  // Fsyn & Fic are in erg/s/eV ; sum is in ergs
	  ssum += (theShocks[i]->Fsyn(m_energy[en],tt*m_dt)+cst.flagIC*theShocks[i]->Fic(m_energy[en],tt*m_dt))*m_de[en]*m_dt;
	}
      }
      theShocks[i]->setSum(ssum);
      //theShocks[i]->Write();
    }
}

/*------------------------------------------------------*/
void GRBSim::ComputeFlux(double time){
  double nshock=theShocks.size();
  double norma;
  double sum;
  double temp;
  m_spectrum.clear();
  
  for (int en=0;en<cst.enstep;en++)
    {
      m_spectrum.push_back(0.0);
    }
  
  for (int j=0;j<nshock;j++)
    {
      sum = theShocks[j]->Sum(); /// erg
      norma = (theShocks[j]->Eint())/(myUni->Area()); //erg/cm^2
      for (int en=0;en<cst.enstep;en++)
	if(sum>0.0)
	  {
	    // temp is in erg/s/eV
	    temp=(theShocks[j]->Fsyn(m_energy[en],time)+
	      cst.flagIC*theShocks[j]->Fic(m_energy[en],time));
	    
	    m_spectrum[en]+=(cst.erg2MeV*1.0e+16)/m_energy[en]*(norma/sum)*temp;
	    // (eV/erg) (1/eV) (erg/cm^2) (1/erg) (erg/s/eV) = 1/cm^2/s/eV
	    // converted in ->photons/s/MeV/m^2
	    
	    //  While...
	    //m_spectrum[en]+=(norma/sum)*temp;
	    // erg/s/eV/cm^2
	  }
    }
}
/*------------------------------------------------------*/
double GRBSim::IFlux(double enmin){
  // gives the integrated flux of energy (eV/s/m^2) for energy > enmin 
  double flux=0.0;
  for (int en=0;en<cst.enstep;en++)
    {
      if(m_energy[en]>=enmin){
	flux +=m_energy[en]*m_spectrum[en]*(m_de[en])*(1.0e-6);
	//eV/s/m^2
      }
    }
  return flux;
}
/*------------------------------------------------------*/
double GRBSim::IRate(double enmin){
  // flux of particles for energy > enmin
  double flux=0.0;
  for (int en=0;en<cst.enstep;en++)
    {
      if(m_energy[en]>=enmin){  
	flux += m_spectrum[en]*(m_de[en])*(1.0e-6);
	//ph/s/m^2
      }
    }
  return flux;  
}
/*------------------------------------------------------*/
double GRBSim::IEnergy(double enmin){
  // Integrated flux of energy for energy > enmin
  // that flows in 1 m_dt .
  double flux=0.0;
  for (int en=0;en<cst.enstep;en++)
    {
      if(m_energy[en]>=enmin){  
	flux +=m_energy[en]*m_spectrum[en]*(m_de[en])*(1.0e-6)*m_dt;
	//eV/m^2
      }
    }
  return flux;  
}
/*------------------------------------------------------*/
float GRBSim::DrawPhotonFromSpectrum(std::vector<double> spctrmVec,double emin)
{
  // Return the energy of a photon using the spectrum.
  // The photon has energy > emin  
  if(spctrmVec.size()) { 
    TH1D *hh1 = new TH1D("hh1","hh1",cst.enstep,log10(cst.enmin),log10(cst.enmax));
    hh1->Reset();
    for (int en=0;en<cst.enstep;en++)
      if (m_energy[en]>=emin) {
	hh1->SetBinContent(en+1,spctrmVec[en]);
      } else {
	hh1->SetBinContent(en+1,0.0);
      }
    double ph = hh1->GetRandom();
    delete hh1;
    return pow(10,ph-9.0); //ph is in log scale... and NRJ in GeV
  } else {
    return 0.0;
  }
}

/*------------------------------------------------------*/
void GRBSim::TotalProperties(int flag){
  double ftottot=0;
  std::vector<double>my_spec;
  TH1D *h1 = new TH1D("h1","h1",cst.enstep,log10(cst.enmin),log10(cst.enmax));
  TCanvas *cc2 = new TCanvas();
  h1->GetXaxis()->SetTitle("Log(Energy [eV])");
  h1->GetXaxis()->SetLabelSize(0.05);
  h1->GetXaxis()->SetTitleOffset(1.5);
  //  h1->GetYaxis()->SetTitle("flux[ph/s/eV/m^2]");  
  h1->GetYaxis()->SetTitle("flux[erg/s/cm^2]");  
  //  h1->SetFillColor(5);

  for (int t=0;t<cst.nstep;t++){
    for (int en=0;en<cst.enstep;en++){
      m_flux[en][t]=0.0;
    }
  }
  
  for (int t=0;t<cst.nstep;t++){
    m_time[t]=t*m_dt;
    // Compute the flux @ time
    ComputeFlux(m_time[t]);
    my_spec=Spectrum(); // is in photons/s/MeV/m^2
    for (int en=0;en<cst.enstep;en++){
      m_flux[en][t]=my_spec[en]*m_energy[en]/(cst.erg2MeV*1.0e+6);  
      // m_flux is in erg/s/MeV/m^2
    }
    if (flag==1) {
      for (int en=0;en<cst.enstep;en++){
	h1->SetBinContent(en+1,m_flux[en][t]*1.0e-4); //erg/cm^2/s/MeV
	//	h1->SetBinContent(en+1,Flux(en)); //ph/m^2/s/eV
      }
      cc2->cd();
      h1->Draw("AL");
      h1->GetXaxis()->SetTitle("Log(Energy[eV])");
      h1->GetXaxis()->SetTitleSize(0.035);
      h1->GetXaxis()->SetLabelSize(0.035);
      h1->GetXaxis()->SetTitleOffset(1.4);
      h1->GetYaxis()->SetTitle("Fv[erg/cm^2/s/MeV]");
      h1->GetYaxis()->SetTitleSize(0.035);
      h1->GetYaxis()->SetLabelSize(0.035);
      h1->GetYaxis()->SetTitleOffset(1.4);   
      
      h1->SetLineWidth(3);
      h1->SetLineColor(4);
      cc2->Update();   
    }
    ftottot+=IEnergy(cst.enmin)/(cst.erg2MeV*1.0e+10)*myUni->Area();
    //erg
    cout<<"Time / tmax = "<<m_time[t]<<"/" <<m_tmax<<" Ftot [erg]= "<<ftottot<<endl;
  }
  cout<<"Fluence [erg/cm^2] = "<<ftottot/myUni->Area()<<endl;
  
}
/*------------------------------------------------------*/
void GRBSim::Plot(){
  
  TCanvas *c1 = new TCanvas("c1","GRB Flux");
  TCanvas *c1b = new TCanvas("c1b","Photons Count");
  TCanvas *c2 = new TCanvas("c2","GRB Light Curve");
  TCanvas *c4 = new TCanvas("c4","Countour Plot");
  TCanvas *c5 = new TCanvas("c5","3D Flux Rapresentation");
  cout<<"xmin = "<< log10(m_energy[0]) << "xmax = "<<m_energy[cst.enstep]<< "ymin = "<<m_time[0] <<"ymax =" << m_time[cst.nstep-1] <<endl;
  
  c1->SetLogx();
  c1->SetLogy();
  c1->SetGridx();
  c1->SetGridy();

  c1b->SetLogx();
  c1b->SetLogy();
  c1b->SetGridx();
  c1b->SetGridy();

  
  TGraph *gfl = new TGraph(cst.enstep);
  TGraph *gph = new TGraph(cst.enstep);

  TGraph *glc = new TGraph(cst.nstep);
  TGraph *gch1 = new TGraph(cst.nstep);
  TGraph *gch2 = new TGraph(cst.nstep);
  TGraph *gch3 = new TGraph(cst.nstep);
  TGraph *gch4 = new TGraph(cst.nstep);
  
  TH2D *h2 = new TH2D("h2","H2",cst.enstep-2,log10(m_energy[1]),log10(m_energy[cst.enstep]),cst.nstep-1,m_time[0],m_time[cst.nstep-1]);
  //TH2D *h2 = new TH2D("h2","H2",cst.enstep-1,2,20,cst.nstep-1,0,1);
  //h1->SetOption("B");
  
  double Fv[cst.enstep];
  double ph_m2_s[cst.enstep];
  double lct[cst.nstep];
  double ch1[cst.nstep];
  double ch2[cst.nstep];
  double ch3[cst.nstep];
  double ch4[cst.nstep];
  
  for (int t=0;t<cst.nstep;t++){
    ch1[t]=0.0;
    ch2[t]=0.0;
    ch3[t]=0.0;
    ch4[t]=0.0;
    lct[t]=0.0;
  }
  double loge[cst.enstep];
  for (int en=0;en<cst.enstep;en++){
    Fv[en]=0.0;
    ph_m2_s[en]=0.0;
    loge[en]=log10(m_energy[en]);    
    for (int t=0;t<cst.nstep;t++) {
      m_flux[en][t]<=1e-25?(m_flux[en][t]=1.0e-25):(m_flux[en][t]);
     
      Fv[en]+=(m_flux[en][t])*/*m_de[en]*/m_dt;  /// is in erg/cm^2/MeV
      ph_m2_s[en]+=(m_flux[en][t])*(cst.erg2MeV*1.0e+6)/m_energy[en]*m_dt;/*m_de[en];*/
      // is in ph/MeV/m^2/s
      lct[t]+=(m_flux[en][t])*(cst.erg2MeV)*m_de[en]/m_energy[en]*1.0e-4; /// is in ph/s/cm^2
      
      if (m_energy[en]<cst.ch1H &&m_energy[en]>=cst.ch1L){
	ch1[t]+=(m_flux[en][t])*(cst.erg2MeV)*m_de[en]/m_energy[en]*1.0e-4;
      }
      if (m_energy[en]<cst.ch2H &&m_energy[en]>=cst.ch2L){
	ch2[t]+=(m_flux[en][t])*(cst.erg2MeV)*m_de[en]/m_energy[en]*1.0e-4;
      }
      if (m_energy[en]<cst.ch3H &&m_energy[en]>=cst.ch3L){
	ch3[t]+=(m_flux[en][t])*(cst.erg2MeV)*m_de[en]/m_energy[en]*1.0e-4;
      }
      if (m_energy[en]<cst.ch4H &&m_energy[en]>=cst.ch4L){
	ch4[t]+=(m_flux[en][t])*(cst.erg2MeV)*m_de[en]/m_energy[en]*1.0e-4;
      }
      h2->Fill(loge[en],m_time[t],m_flux[en][t]);
      //      cout<<log10(m_energy[en])<<m_time[t]<<log10(m_flux[en][t])<<endl;
    }
  }
  
  gfl=new TGraph(cst.enstep,m_energy,Fv);
  gph=new TGraph(cst.enstep,m_energy,ph_m2_s);

  glc=new TGraph(cst.nstep,m_time,lct);
  
  gch1=new TGraph(cst.nstep,m_time,ch1);
  gch2=new TGraph(cst.nstep,m_time,ch2);
  gch3=new TGraph(cst.nstep,m_time,ch3);
  gch4=new TGraph(cst.nstep,m_time,ch4);
  
  gfl->SetLineWidth(2);
  glc->SetLineWidth(2);
  gfl->SetLineColor(4);

  gch1->SetLineColor(5);
  gch2->SetLineColor(3);
  gch3->SetLineColor(4);
  gch4->SetLineColor(2);
  
  
  c1->cd();
  gfl->Draw("ALP");
  gfl->GetXaxis()->SetTitle("Energy [eV]");
  gfl->GetXaxis()->SetTitleSize(0.035);
  gfl->GetXaxis()->SetLabelSize(0.035);
  gfl->GetXaxis()->SetTitleOffset(1.4);
  gfl->GetYaxis()->SetTitle("Fv[erg/cm^2/MeV]");
  gfl->GetYaxis()->SetTitleSize(0.035);
  gfl->GetYaxis()->SetLabelSize(0.035);
  gfl->GetYaxis()->SetTitleOffset(1.4);
  gfl->Draw("ALP");
  
  c1b->cd();
  gph->Draw("ALP");
  gph->GetXaxis()->SetTitle("Log(Energy [eV])");
  gph->GetXaxis()->SetTitleSize(0.035);
  gph->GetXaxis()->SetLabelSize(0.035);
  gph->GetXaxis()->SetTitleOffset(1.4);
  gph->GetYaxis()->SetTitle("Particle/cm^2/s/MeV]");
  gph->GetYaxis()->SetTitleSize(0.035);
  gph->GetYaxis()->SetLabelSize(0.035);
  gph->GetYaxis()->SetTitleOffset(1.4);
  gph->Draw("ALP");
  
  c1b->Update();

  c2->cd();
  TLegend *leg = new TLegend(0.8,0.8,0.98,0.98);
  leg->AddEntry(gch1,"ch1","l");
  leg->AddEntry(gch2,"ch2","l");
  leg->AddEntry(gch3,"ch3","l");
  leg->AddEntry(gch4,"ch4","l");
  leg->AddEntry(glc,"Sum","l");
  
  glc->Draw("ALP");
  glc->GetXaxis()->SetTitle("Time [sec]");
  glc->GetXaxis()->SetTitleSize(0.035);
  glc->GetXaxis()->SetLabelSize(0.035);
  glc->GetXaxis()->SetTitleOffset(1.4);
  glc->GetYaxis()->SetTitle("flux[ph/cm^2/sec]");
  glc->GetYaxis()->SetTitleSize(0.035);
  glc->GetYaxis()->SetLabelSize(0.035);
  glc->GetYaxis()->SetTitleOffset(1.4);
  glc->Draw("ALP");
  gch1->Draw("LP");
  gch2->Draw("LP");
  gch3->Draw("LP");
  gch4->Draw("LP");
  leg->Draw();
  c2->Update();
  
  c4->cd(); 
  h2->Draw("CONT");
  c4->Update();
 
  c5->cd();
  h2->Draw("surf");
  h2->GetXaxis()->SetTitle("Time");
  //  h2->GetYaxis()->SetTitleSize(0.05);
  h2->GetXaxis()->SetLabelSize(0.05);
  h2->GetXaxis()->SetTitleOffset(1.5);
  h2->GetYaxis()->SetTitle("Log(Energy [eV])");
  //h2->GetYaxis()->SetTitleSize(0.035);
  h2->GetYaxis()->SetLabelSize(0.05);
  h2->GetYaxis()->SetTitleOffset(1.5);
  h2->GetZaxis()->SetTitle("vFv[erg/s/m^2/MeV]");
  // h2->GetYaxis()->SetTitleSize(0.035);
  h2->GetZaxis()->SetLabelSize(0.05);
  h2->GetZaxis()->SetTitleOffset(1.3);
  h2->Draw("surf");
  c5->Update();
};













