#include <iomanip>
#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TStyle.h"

//gROOT->Reset();
TFile myFile;
TFile myOutFile;
TTree *myTree;
TString TID;

const int NUM_ROC = 2;
const int NUM_LAYER = 18;
const int NUM_VIEW = 2;

void help();


void AddTray(int l, int v)
{
  std::cout<<TrayId(l,v)<<" add to the stack of tree data "<<std::endl;
  myTree->AddFriend(TrayId(l,v),&myFile);
}

void SelectTray(int l, int v)
{
  std::cout<<TrayId(l,v)<<" Selected for plotting "<<std::endl;
}

void Init(char *filename="MyRootFile.root", char* f="MyAnalysisFile.root")
{
  
  gStyle->SetCanvasColor(10);
  gStyle->SetFillColor(10);
  
  TFile myFile(filename);
  
  myTree = (TTree*) myFile.Get("Header");

  AddTray(1,0);
  SelectTray(1,0);

  TFile myOutFile(f, "recreate");
  help();
}


void help()
{
  std::cout<<" ******************** Analysis Program for Tray data ********************"<<std::endl;
  std::cout<<" ---------- DIGI ANALYSIS:     "<<std::endl;
  std::cout<<" Init(input filename, outFilename) : Initialize the analysis. "<<endl;
  std::cout<<" AddTray(l,d) : Add a tray to the stack of the trays to analyze. 
                              By default Tray(1,0) is added                 "<<std::endl;
  std::cout<<" SelectTray(l,d): The Tray l,d is selected for plotting       "<<std::endl;
  std::cout<<" TrayId(l,d) "<<std::endl;
  std::cout<<" PlotToT(TString myCuts='')  "<<std::endl;
  std::cout<<" PlotTotTot(TString myCuts='')"<<std::endl;
  std::cout<<" PlotHitMap(TString myCuts='') "<<std::endl;
  std::cout<<" PlotNumHitsLayer(TString myCuts='') "<<std::endl;
  std::cout<<" FindEvents(TString myCuts='') "<<std::endl;
  std::cout<<""<<std::endl;
  std::cout<<" ---------- RECON ANALYSIS (only if available):     "<<std::endl;
  std::cout<<" PlotVtx(TString myCuts='')"<<std::endl;
  std::cout<<" PlotClusHit(TString myCuts='')"<<std::endl;
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


void PlotToT_0(TString myCuts="")
{
  
  TString var="ToT0";
  var.Prepend(TID);

  gDirectory->Delete("ToT0");
  TH1D *ToT0 = new TH1D("ToT0",var,255,0,255);
  var+=">>ToT0";
  myTree->Draw(var,myCuts);
}

void PlotToT_1(TString myCuts="")
{
  TString var="ToT1";
  var.Prepend(TID);

  gDirectory->Delete("ToT1");
  TH1D *ToT1 = new TH1D("ToT1",var,255,0,255);
  var+=">>ToT1";
  myTree->Draw(var,myCuts);
}

void PlotToT(TString myCuts="")
{
  gDirectory->Delete("ToTC");
  TCanvas *ToTC = new TCanvas("ToTC","ToTC",800,400);
  ToTC->Divide(2,1);
  ToTC->cd(1);
  PlotToT_0(myCuts);
  ToTC->cd(2);
  PlotToT_1(myCuts);
}

void PlotHitMap(TString myCuts="")
{
  gDirectory->Delete("HitCan");
  TCanvas *HitCan = new TCanvas("HitCan","Hit Map",1200,500);
  TString var="TkrHits";
  var.Prepend(TID);
  gDirectory->Delete("HitMap");
  TH1D *HitMap = new TH1D("HitMap",var,1536,0,1536);
  HitMap->SetLineColor(4);
  var+=">>HitMap";
  myTree->Draw(var,myCuts);
  std::cout<<" Number of Hits: "<<HitMap->GetEntries()<<std::endl;
}


void ProjectHitMap(TString myCuts="")
{
  gDirectory->Delete("prjHitCan");
  TCanvas *prjHitCan = new TCanvas("prjHitCan","Projected Hit Map",500,500);
  TString var="TkrHits";
  var.Prepend(TID);
  gDirectory->Delete("HitMap");
  TH1D *HitMap = new TH1D("HitMap",var,1536,0,1536);

  var+=">>HitMap";
  myTree->Draw(var,myCuts);

  int max_hit=HitMap->GetMaximum();
  //  int min_hit=HitMap->GetMinimum();
  
  TH1D *PrjHitMap = new TH1D("PrjHitMap","Projected Hit Map",100,0,max_hit);
  
  for(int strip=0;strip<1536;)
    {
      PrjHitMap->Fill(HitMap->GetBinContent(strip++));      
    }


  PrjHitMap->SetFillColor(5);
  PrjHitMap->Draw();

}





void PlotNumHitsLayer(TString myCuts="")
{
  gDirectory->Delete("NumHitsCan");
  TCanvas *NumHitsCan = new TCanvas("NumHitsCan","Num Hits",500,500);
  TString var="TkrNumHits";
  var.Prepend(TID);
  gDirectory->Delete("TkrNumHits");
  //  TH1D *TkrNumHits = new TH1D("TkrNumHits",var,1536,0,1536);
  var+=">>TkrNumHits";
  myTree->Draw(var,myCuts);
  TkrNumHits->SetFillColor(3);
  std::cout<<" Number of Events: "<<myTree->GetEntries()<<std::endl;
}

void FindEvents(TString myCuts="")
{
  myTree->Scan("EventId",myCuts);
}

void PlotTotTot(TString myCuts="")
{ 
  gDirectory->Delete("ToTToTCan");
  TCanvas *ToTToTCan = new TCanvas("ToTToTCan","ToT vs ToT",500,500);
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

void PlotAllDigi(TString myCuts="")
{
  PlotToT(myCuts);
  PlotTotTot(myCuts);
  PlotHitMap(myCuts);
  PlotNumHitsLayer(myCuts);
}
//////////////////////////////////////////////////
// RECON
//////////////////////////////////////////////////

void PlotVtx(TString myCuts="")
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

  TCanvas *CanVtx = new TCanvas("CanVtx","CanVtx",1000,300);
  CanVtx->Divide(3,1);
  CanVtx->cd(1);
  myTree->Draw("TkrVtxX",myCuts);
  CanVtx->cd(2);
  myTree->Draw("TkrVtxY",myCuts);
  CanVtx->cd(3);
  myTree->Draw("TkrVtxZ",myCuts);
  
}

void PlotClusHit(TString myCuts="")
{
  TString var0="TkrNumClus";
  TString var1="TkrTotalNumHits";
  TString var=var0;
  var+=":";
  var+=var1;
  gDirectory->Delete("NumClus_NumHits");
  TH2D *Histo = new TH2D("NumClus_NumHits","Clusters vs Hits",30,0,15,30,0,15);
  Histo->SetXTitle("Number of Hits");
  Histo->SetYTitle("Number of Clusters");
  
  var+=">>NumClus_NumHits";
  myTree->Draw(var,myCuts);
  Histo->Draw("box");
  /*
    TPaveStats *st = (TPaveStats*) Histo->GetListOfFunctions()->FindObject("stats");
    st->SetOptStat(11);
  */
  gStyle->SetOptStat(11);
  Histo->Draw("box");

}

