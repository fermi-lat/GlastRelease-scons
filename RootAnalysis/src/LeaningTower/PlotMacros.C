#include <iomanip>
#include <iostream>
#include "TF1.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TTree.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TString.h"
#include "TStyle.h"

TFile *myFile;
TTree *myTree;
TString TID;
////////////////// Functions definitions ////////////////////////////////
void help();
void AddTray(int l, int v);
void SelectTray(int l, int v);
void Init(char *filename="MyRootFile.root");
TString TrayId(long int l, long int v);
void PlotToT_0(TString myCuts="");
void PlotToT_1(TString myCuts="");
void PlotToT(TString myCuts="");
void PlotHitMap(TString myCuts="");
void ProjectHitMap(TString myCuts="");
void PlotNumHitsLayer(TString myCuts="");
void FindEvents(TString myCuts="");
void PlotTotTot(TString myCuts="");
void PlotAllDigi(TString myCuts="");
void PlotVtx(TString myCuts="");
void PlotClusHit(TString myCuts="");

////////////////// Implementation: ////////////////////
void AddTray(int l, int v)
{
  std::cout<<TrayId(l,v)<<" add to the stack of tree data "<<std::endl;
  myTree->AddFriend(TrayId(l,v),myFile);
}

void SelectTray(int l, int v)
{
  AddTray(l,v);
  std::cout<<" & Selected for plotting "<<std::endl;
}

void Init(char *filename)
{
  
  gStyle->SetCanvasColor(10);
  gStyle->SetFillColor(10);
  
  TFile *myFile = new TFile(filename);
  
  myTree = (TTree*) myFile->Get("Header");

  AddTray(1,0);
  SelectTray(1,0);

  //  TFile myOutFile(f, "recreate");
  help();
}


void help()
{
  std::cout<<" ******************** Analysis Program for Tray(s) data ********************"<<std::endl;
  std::cout<<" ---------- DIGI ANALYSIS:     "<<std::endl;
  std::cout<<" Init(input filename, outFilename) : Initialize the analysis. "<<std::endl;
  std::cout<<" AddTray(l,d) : Add a tray to the stack of the trays to analyze.  "<<std::endl;
  std::cout<<"               By default Tray(1,0) is added                     "<<std::endl;
  std::cout<<" SelectTray(l,d): The Tray l,d is selected for plotting       "<<std::endl;
  std::cout<<"                  By default Tray(1,0) is added               "<<std::endl;
  std::cout<<" ------------------------------ Plot Macros: --------------------"<<std::endl;
  std::cout<<" PlotToT(TString myCuts=\"\")  "<<std::endl;
  std::cout<<" PlotTotTot(TString myCuts=\"\")"<<std::endl;
  std::cout<<" PlotHitMap(TString myCuts=\"\") "<<std::endl;
  std::cout<<" ProjectedHitMap(TString myCuts=\"\") "<<std::endl;
  std::cout<<" PlotNumHitsLayer(TString myCuts=\"\") "<<std::endl;
  std::cout<<" FindEvents(TString myCuts=\"\") "<<std::endl;
  std::cout<<" --------------------"<<std::endl;
  std::cout<<" \"myCuts\" can be of the form:  "<<std::endl;
  std::cout<<" \"ToT0>0\" ,\"ToT0>0 && TkrNumHits >10\", \"EventId == 265\"..."<<std::endl;
  std::cout<<" If more than one tray has been added to the stack,       "<<std::endl;
  std::cout<<" than it is possible to specify the tray used for cutting "<<std::endl;
  std::cout<<" (otherwise the selected tray for plotting is used)       "<<std::endl;
  std::cout<<""<<std::endl;
  std::cout<<" For Example:"<<std::endl;
  std::cout<<" PlotToT(\"ToT0>0\")                 Plot the ToT for the Tray2 (defaul) if its ToT is >0"<<std::endl;
  std::cout<<" AddTray(1,1)                        (this add the tray 1,1 to the the stack as Tray3 (do it only once!) "<<std::endl;
  std::cout<<" PlotToT(\"Tray3.ToT0>0\")           Plot the ToT for the Tray2 if the ToT of the Tray3 is is >0"<<std::endl;
  std::cout<<" SelectTray(1,1)                     This select the tray 3..."<<std::endl;
  std::cout<<" PlotToT(\"Tray3.ToT0>0\")           Plot the ToT for the Tray3 if the ToT of the Tray3 is is >0 (same as PlotToT())"<<std::endl;
  std::cout<<" PlotToT(\"Tray2.ToT0>0\")           Plot the ToT for the Tray3 if the ToT of the Tray2 is is >0 (same as PlotToT())"<<std::endl;
  std::cout<<" Note that PlotToT(\"TkrHits[]==5\") Plot the ToTs for the events which hit the strip 5 !!"<<std::endl; 
  std::cout<<" "<<std::endl;
  std::cout<<" The available variables are:"<<std::endl;
  std::cout<<" EventId, RunId, TkrTotalNumHits (Global variables)"<<std::endl;
  std::cout<<" TkrNumHits,  ToT0, ToT1, TkrHits[]  (These are for each trays"<<std::endl;
  std::cout<<" NOTE: Is alwais possible accessing directly to the tree : \"myTree->Draw(\"RunId\"), or myTree->Draw(\"Tray2.ToT0\")"<<std::endl;
  std::cout<<" ---------- RECON ANALYSIS (only if available):     "<<std::endl;
  std::cout<<" PlotVtx(TString myCuts=\"\")"<<std::endl;
  std::cout<<" PlotClusHit(TString myCuts=\"\")"<<std::endl;
  std::cout<<" PlotAllDigi()"<<std::endl;
  
}

//////////////////////////////////////////////////
//              1 Tray Analysis                 //
//////////////////////////////////////////////////
TString TrayId(long int l, long int v)
{
  //  std::cout<<"---------- Layer = "<<l<<" ; view = "<<v<<std::endl;
  int Plane = 2*l+v;
  TString ts="Tray";//"[";
  //  TString *trayid;
  //ts +=l ;
  //ts+="][";
  //ts+=v;
  //ts+="]";
  ts+=Plane;
  TID=ts+".";
  return ts;
}


void PlotToT_0(TString myCuts)
{
  
  TString var="ToT0";
  var.Prepend(TID);

  gDirectory->Delete("ToT0");
  TH1D *ToT0;
  ToT0 = new TH1D("ToT0",var,255,0,255);
  var+=">>ToT0";
  myTree->Draw(var,myCuts);
}

void PlotToT_1(TString myCuts)
{
  TString var="ToT1";
  var.Prepend(TID);

  gDirectory->Delete("ToT1");
  TH1D *ToT1;
  ToT1 = new TH1D("ToT1",var,255,0,255);
  var+=">>ToT1";
  myTree->Draw(var,myCuts);
}

void PlotToT(TString myCuts)
{
  gDirectory->Delete("ToTC");
  TCanvas *ToTC;
  ToTC = new TCanvas("ToTC","ToTC",800,400);
  ToTC->Divide(2,1);
  ToTC->cd(1);
  PlotToT_0(myCuts);
  ToTC->cd(2);
  PlotToT_1(myCuts);
}

void PlotHitMap(TString myCuts)
{
  gDirectory->Delete("HitCan");
  TCanvas *HitCan;
  HitCan = new TCanvas("HitCan","Hit Map",1200,500);
  TString var="TkrHits";
  var.Prepend(TID);
  gDirectory->Delete("HitMap");
  TH1D *HitMap;
  HitMap = new TH1D("HitMap",var,1535,-0.5,1534.5);
  HitMap->SetLineColor(4);
  var+=">>HitMap";
  myTree->Draw(var,myCuts);
  std::cout<<" Number of Hits: "<<HitMap->GetEntries()<<std::endl;
}


void ProjectHitMap(TString myCuts)
{
  gDirectory->Delete("prjHitCan");
  TCanvas *prjHitCan;
  prjHitCan = new TCanvas("prjHitCan","Projected Hit Map",500,500);
  TString var="TkrHits";
  var.Prepend(TID);
  gDirectory->Delete("HitMap1");
  TH1D *HitMap1;
  HitMap1 = new TH1D("HitMap1",var,1535,-0.5,1534.5);
  var+=">>HitMap1";
  myTree->Draw(var,myCuts);

  int max_hit=(int) HitMap1->GetMaximum();
  TH1D *PrjHitMap;
  PrjHitMap = new TH1D("PrjHitMap","Projected Hit Map",100,0,max_hit);
  
  for(int strip=1;strip<=1535;strip++)
    {
      PrjHitMap->Fill(HitMap1->GetBinContent(strip));      
    }
  
  PrjHitMap->SetFillColor(5);
  int Integral = (int) PrjHitMap->GetSum();
  if (Integral!=0)  PrjHitMap->Scale(1./Integral);

  PrjHitMap->Draw();
  PrjHitMap->Fit("gaus","qsame");

  TF1 *Gaussian = PrjHitMap->GetFunction("gaus");
  
  double sigma = Gaussian->GetParameter(2);
  double mean  = Gaussian->GetParameter(1);
  std::cout<<" Mean : "<<mean<<", Sigma "<<sigma<<std::endl;
  
  for(int strip=1;strip<=1535;strip++)
    {
      int numHits = (int) HitMap1->GetBinContent(strip);
      if(numHits < mean - 3.0*sigma)
	std::cout<<" WARNING!! Strip Number "<<strip-1<<" is below 3 sigma from the mean!! Counts : "<<numHits<<std::endl;
      
      if(numHits > mean + 3.0*sigma)
	std::cout<<" WARNING!! Strip Number "<<strip-1<<" is above 3 sigma from the mean!! Counts : "<<numHits<<std::endl;
    }
  
  
}

void PlotNumHitsLayer(TString myCuts)
{
  gDirectory->Delete("NumHitsCan");
  TCanvas *NumHitsCan;
  NumHitsCan = new TCanvas("NumHitsCan","Num Hits",500,500);
  TString var="TkrNumHits";
  var.Prepend(TID);
  gDirectory->Delete("TkrNumHits");
  var+=">>TkrNumHits";
  myTree->Draw(var,myCuts);
  TH1D *TkrNumHits = (TH1D*) gDirectory->Get("TkrNumHits");
  TkrNumHits->SetFillColor(3);
  std::cout<<" Number of Events: "<<myTree->GetEntries()<<std::endl;
}

void FindEvents(TString myCuts)
{
  myTree->Scan("EventId",myCuts);
}

void PlotTotTot(TString myCuts)
{ 
  gDirectory->Delete("ToTToTCan");
  TCanvas *ToTToTCan;
  ToTToTCan = new TCanvas("ToTToTCan","ToT vs ToT",500,500);
  TString var0="ToT0";
  TString var1="ToT1";
  var0.Prepend(TID);
  var1.Prepend(TID);
  TString var=var0;
  var+=":";
  var+=var1;
  
  //  var+=">>tot";
  myTree->Draw(var,myCuts);
  TGraph *graph = (TGraph*) gPad->GetPrimitive("Graph");
  graph->SetMarkerStyle(25);
  graph->SetMarkerSize(0.2);
  graph->SetMarkerColor(2);
  graph->Draw("P");
  
}

void PlotAllDigi(TString myCuts)
{
  PlotToT(myCuts);
  PlotTotTot(myCuts);
  PlotHitMap(myCuts);
  PlotNumHitsLayer(myCuts);
  ProjectHitMap(myCuts);
}

void DumpReport()
{
  TCanvas ReportCanvas("ReportCanvas");
  ReportCanvas.Print("Report.ps[");

  



  TCanvas *SelectedCanvas;
  SelectedCanvas = (TCanvas*) gROOT->FindObject("ToTC");
  if(SelectedCanvas) SelectedCanvas->Print("Report.ps");
  SelectedCanvas = (TCanvas*) gROOT->FindObject("HitCan");
  if(SelectedCanvas) SelectedCanvas->Print("Report.ps");
  SelectedCanvas = (TCanvas*) gROOT->FindObject("prjHitCan");
  if(SelectedCanvas) SelectedCanvas->Print("Report.ps");
  SelectedCanvas = (TCanvas*) gROOT->FindObject("NumHitsCan");
  if(SelectedCanvas) SelectedCanvas->Print("Report.ps");
  SelectedCanvas = (TCanvas*) gROOT->FindObject("ToTToTCan");
  if(SelectedCanvas) SelectedCanvas->Print("Report.ps");
  
  ReportCanvas.Print("Report.ps]");
}
//////////////////////////////////////////////////
// RECON
//////////////////////////////////////////////////

void PlotVtx(TString myCuts)
{
  //  TString vtxX="TkrVtxX";
  
  TString NumVtx="TkrNumVtx > 0";

  if(myCuts!="")
    {
      myCuts+=" & ";
      myCuts+=NumVtx;
    }
  else
    myCuts=NumVtx;

  TCanvas *CanVtx;
  CanVtx = new TCanvas("CanVtx","CanVtx",1000,300);
  CanVtx->Divide(3,1);
  CanVtx->cd(1);
  myTree->Draw("TkrVtxX",myCuts);
  CanVtx->cd(2);
  myTree->Draw("TkrVtxY",myCuts);
  CanVtx->cd(3);
  myTree->Draw("TkrVtxZ",myCuts);
  
}

void PlotClusHit(TString myCuts)
{
  TString var0="TkrNumClus";
  TString var1="TkrTotalNumHits";
  TString var=var0;
  var+=":";
  var+=var1;
  gDirectory->Delete("NumClus_NumHits");
  TH2D *NumClus_NumHits;
  NumClus_NumHits = new TH2D("NumClus_NumHits","Clusters vs Hits",30,0,15,30,0,15);
  NumClus_NumHits->SetXTitle("Number of Hits");
  NumClus_NumHits->SetYTitle("Number of Clusters");
  
  var+=">>NumClus_NumHits";
  myTree->Draw(var,myCuts);
  NumClus_NumHits->Draw("box");
  /*
    TPaveStats *st = (TPaveStats*) NumClus_NumHits->GetListOfFunctions()->FindObject("stats");
    st->SetOptStat(11);
  */
  gStyle->SetOptStat(11);
  NumClus_NumHits->Draw("box");

}

