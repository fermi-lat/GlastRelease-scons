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
#include "TPaveText.h"

TFile *myFile;
TTree *myTree;
TString TID;
TString LayerName;
long NumberOfEvents;
////////////////// Functions definitions ////////////////////////////////
void help();
void AddRecon();
void AddLayer(TString VL);
void Initialize(TString filename="MyRootFile.root");
TString LayerId(TString VL);

////////// HISTOGRAMS //////////
TH1D *ToT_0(TString myCuts="");
TH1D *ToT_1(TString myCuts="");
TH1D *HitMap(TString myCuts="");
TH1D *ProjectHitMap(TString myCuts="");
TH1D *NumHitsLayer(TString myCuts="");
TGraph *Correlation(TString var, TString VL1, TString VL2, TString myCuts="");
TGraph *Graph(TString VL1,  TString varA, TString varB, TString myCuts="");
TGraph *TimeHistory(TString VL1);

////////// CANVASES ////////////
void PlotToT(TString myCuts="");
void PlotHitMap(TString myCuts="");
void PlotProjectHitMap(TString myCuts="");
void PlotNumHitsLayer(TString myCuts="");

void PlotCorrelation(TString var, TString VL1, TString VL2, TString myCuts="");
void PlotGraph(TString VL1,  TString varA, TString varB, TString myCuts="");
void PlotHitsVsTime(TString VL1);
//--------------------------------------------------//
long FindEvents(TString myCuts="");
long FindBigEvents(int cut);
void  PlotVtx(TString myCuts="");
TH2D *ClusHit(TString myCuts="");
void PlotClusHit(TString myCuts="");
void PlotRecon(TString myCuts="");
//--------------------------------------------------//
void PlotAllDigi(TString LV="X0",TString myCuts="");
TCanvas *Report(TString LV="X0",TString myCuts="");
void ReportAll(TString="", TString="gif");
void DumpReport(TString LV="X0",char *filename="MyRootFile.root");
////////////////// Implementation: ////////////////////
void AddLayer(TString VL)
{
  std::cout<<LayerId(VL)<<" add to the stack of tree data "<<std::endl;
  myTree->AddFriend(LayerId(VL),myFile);
}

void AddRecon()
{
  myTree->AddFriend("Recon",myFile);
  std::cout<<" ... if everithing is ok, now you have also the following variables:"<<std::endl;
  std::cout<<" Phi, Theta, ThetaXZ, ThetaYZ, TkrClusPlane, TkrClusView"<<std::endl;
  std::cout<<"  TkrClusX, TkrClusY, TkrClusZ"<<std::endl;
  std::cout<<"  TkrVtxX, TkrVtxY, TkrVtxZ, TktTheta"<<std::endl;
  std::cout<<" TkrNumClus, TkrNumVtx, TkrNumTraks, TkrTrk1Cluster"<<std::endl;
}

void Initialize(TString filename)
{
  
  gStyle->SetCanvasColor(10);
  gStyle->SetFillColor(10);
  
  myFile = new TFile(filename);
  
  myTree = (TTree*) myFile->Get("Header");
  NumberOfEvents = (long) myTree->GetEntries();
  //  myTree->AddFriend("Recon",myFile);
  //  AddLayer("X0");
  //  TFile myOutFile(f, "recreate");

  help();
}


void help()
{
  std::cout<<" ******************** Analysis Program for Layer(s) data ********************"<<std::endl;
  std::cout<<" ---------- DIGI ANALYSIS:     "<<std::endl;
  std::cout<<" 1) Initialize(input filename, outFilename) : Initialize the analysis. "<<std::endl;
  std::cout<<" 1.5) AddRecon()                            : Add the recon tree.  "<<std::endl;
  std::cout<<" 2) AddLayer(TString LV)   : Add a layer to the stack of the layers to analyze and select it for plotting.  "<<std::endl;
  std::cout<<"                    example: AddLayer(\"X1\") Add the LayerX1..."<<std::endl;
  std::cout<<" 3) ------------------------------ Plot Macros: --------------------"<<std::endl;
  std::cout<<" PlotToT(TString myCuts=\"\")  "<<std::endl;
  std::cout<<" PlotHitMap(TString myCuts=\"\") "<<std::endl;
  std::cout<<" ProjectedHitMap(TString myCuts=\"\") "<<std::endl;
  std::cout<<" PlotNumHitsLayer(TString myCuts=\"\") "<<std::endl;
  std::cout<<" FindEvents(TString myCuts=\"\") "<<std::endl;
  std::cout<<" PlotCorrelation(TString var, TString VL1, TString VL2, TString myCuts=\"\");            "<<std::endl;  
  std::cout<<"                    Plot a scatter plot of var of the layer LV1 vs var of the Layer LV2  "<<std::endl;  
  std::cout<<" PlotGraph(TString VL1,  TString varA, TString varB, TString myCuts=\"\")                "<<std::endl;  
  std::cout<<"                    Plot a scatter plot of varA of the layer VL1 vs varB of the Layer VL1"<<std::endl;  
  std::cout<<" PlotHitsVsTime(TString VL1);                                               "<<std::endl;                       
  std::cout<<"                    Plot Hits of the layer VL1 vs EventId                                 "<<std::endl;                          
  std::cout<<" --------------------"<<std::endl;
  std::cout<<" \"myCuts\" can be of the form:  "<<std::endl;
  std::cout<<" \"ToT0>0\" ,\"ToT0>0 && TkrNumHits >10\", \"EventId == 265\"..."<<std::endl;
  std::cout<<" If more than one tray has been added to the stack,       "<<std::endl;
  std::cout<<" than it is possible to specify the tray used for cutting "<<std::endl;
  std::cout<<" (otherwise the selected tray for plotting is used)       "<<std::endl;
  std::cout<<""<<std::endl;
  std::cout<<" READ THIS EXAMPLE OF ANALYSIS:"<<std::endl;
  std::cout<<" Initialize(TString filename=\"MyRootFile.root\") Load the file (or reset the analysis)           "<<std::endl; 
  std::cout<<" AddLayer(\"X0\")                    Select Layer  X0 "<<std::endl;
  std::cout<<" PlotToT(\"ToT0>0\")                  Plot the ToT for the LayerX0 if its ToT0 is >0"<<std::endl;
  std::cout<<" AddLayer(\"X1\")                    Add the Layer X1 to the stack of layers "<<std::endl;
  std::cout<<" PlotToT(\"ToT0>0\")                  Plot the ToT for the LayerX1 if the ToT0 of the Layer X0 is >0"<<std::endl;
  std::cout<<" PlotToT(\"LayerX0.ToT0>0\")          Plot the ToT for the LayerX1 if the ToT0 of the Layer X0 is >0"<<std::endl;
  std::cout<<" PlotToT(\"LayerX1.ToT0>0\")          Plot the ToT for the LayerX1 if its ToT0 is >0"<<std::endl;
  std::cout<<" AddLayer(\"X0\")                    Add (select) the Layer X0"<<std::endl;
  std::cout<<" PlotToT(\"LayerX1.ToT0>0\")          Plot the ToT for the LayerX0 if the ToT of the LayerX1 is is >0"<<std::endl;
  std::cout<<" Note that PlotToT(\"TkrHits[]==5\") Plot the ToTs for the events which hit the strip 5 !!"<<std::endl; 
  std::cout<<" --------------------------------------------------    "<<std::endl; 
  std::cout<<" Initialize(TString filename=\"MyRootFile.root\") Load the file (or reset the analysis)           "<<std::endl; 
  std::cout<<" AddLayer(\"X0\")                    Select Layer  X0 "<<std::endl;
  std::cout<<" PlotToT(\"ToT0>0\")                  Plot the ToT for the LayerX0 if its ToT0 is >0"<<std::endl;
  std::cout<<" Initialize(TString filename=\"MyRootFile.root\") Load the file (This reset the analysis!!)        "<<std::endl; 
  std::cout<<" AddLayer(\"X1\")                    Add the Layer X1 to the stack of layers "<<std::endl;
  std::cout<<" PlotToT(\"ToT0>0\")                  Plot the ToT for the LayerX1 its ToT0 is >0 (different from above)"<<std::endl;
  std::cout<<" "<<std::endl;
  std::cout<<" ------------------------------ The available variables are:------------------------------"<<std::endl;
  std::cout<<" EventId,  RunId,  TkrTotalNumHits (Global variables)"<<std::endl;
  std::cout<<" TkrNumHits,  ToT0, ToT1, TkrHits[]  (These are for each layer)"<<std::endl;
  std::cout<<" NOTE: Is alwais possible accessing directly to the tree : \"myTree->Draw(\"RunId\"), or myTree->Draw(\"LayerX1.ToT0\")"<<std::endl;
  std::cout<<" -----------------------------------------------------------------------------------------"<<std::endl;  
  std::cout<<" ---------- RECON ANALYSIS (only if available):     "<<std::endl;
  std::cout<<" PlotVtx(TString myCuts=\"\")"<<std::endl;
  std::cout<<" PlotClusHit(TString myCuts=\"\")"<<std::endl;
  std::cout<<" PlotRecon(TString myCuts=\"\")"<<std::endl;
  std::cout<<" ====================    Macros    ===================="<<std::endl;
  std::cout<<" 1) Initialize(TString filename=\"MyRootFile.root\") Load the file"<<std::endl; 
  std::cout<<" 2) AddLayer(\"X0\") Add the layer 0X to the stak of layers... ready to be analyzed"<<std::endl; 
  std::cout<<" 3) PlotAllDigi(TString myCuts=\"\") Execute all the Plot macros for digi root file..."<<std::endl;
  std::cout<<" 4) Report(TString myCuts=\"\") Makes a report page with some plots and infos. "<<std::endl;
  std::cout<<" 5) DumpReport(char *filename=\"MyRootFile.root\")   Execute different Report() with different set of cuts (see code)"<<std::endl;
  std::cout<<"                                                     It saves all the plot in a ps file (Report_[layer].ps)"<<std::endl;
  
}

//////////////////////////////////////////////////
//              1 Layer Analysis                 //
//////////////////////////////////////////////////
TString LayerId(TString VL)
{
  LayerName="Layer";
  LayerName+=VL;
  TID=LayerName+".";
  return LayerName;
}


TH1D *ToT_0(TString myCuts)
{
  TString var="ToT0";
  var.Prepend(TID);
  TString title=var;
  if(myCuts!="")
    {
      title+=" {";
      title+=myCuts;
      title+="}";
    }
  //  gDirectory->Delete("ToT0");
  TH1D *ToT0;
  ToT0 = new TH1D("ToT0",title,255,0,255);
  var+=">>ToT0";
  myTree->Draw(var,myCuts);
  return ToT0;
}

TH1D *ToT_1(TString myCuts)
{
  TString var="ToT1";
  var.Prepend(TID);
  TString title=var;
  if(myCuts!="")
    {
      title+=" {";
      title+=myCuts;
      title+="}";
    }
//  gDirectory->Delete("ToT1");
  TH1D *ToT1;
  ToT1 = new TH1D("ToT1",title,255,0,255);
  var+=">>ToT1";
  myTree->Draw(var,myCuts);
  return ToT1;
}

void PlotToT(TString myCuts)
{
  gDirectory->Delete("ToTC");
  TCanvas *ToTC;
  ToTC = new TCanvas("ToTC","ToTC",800,400);
  ToTC->Divide(2,1);
  ToTC->cd(1);
  ToT_0(myCuts);
  ToTC->cd(2);
  ToT_1(myCuts);
}

TH1D *HitMap(TString myCuts)
{
  TString var="TkrHits";
  var.Prepend(TID);
  TString title=var;
  if(myCuts!="")
    {
      title+=" {";
      title+=myCuts;
      title+="}";
    }
  //  gDirectory->Delete("HitMap");
  TH1D *HitMap;
  HitMap = new TH1D("HitMap",title,1536,-0.5,1535.6);
  HitMap->SetLineColor(4);
  var+=">>HitMap";
  myTree->Draw(var,myCuts);
  std::cout<<" Number of Hits: "<<HitMap->GetEntries()<<std::endl;
  return (TH1D*) HitMap;
}


void PlotHitMap(TString myCuts)
{
  gDirectory->Delete("HitCan");
  TCanvas *HitCan;
  HitCan = new TCanvas("HitCan","Hit Map",1200,500);
  HitMap(myCuts);
}


TH1D *ProjectHitMap(TString myCuts)
{
  gDirectory->Delete("myCanvas");
  TCanvas *myCanvas;
  myCanvas = new TCanvas("myCanvas","myCanvas",100,50);

  gDirectory->Delete("HitMap1");
  
  TString var="TkrHits";
  var.Prepend(TID);
  TString title=var;
  if(myCuts!="")
    {
      title+=" {";
      title+=myCuts;
      title+="}";
    }
  gDirectory->Delete("HitMap1");
  TH1D *HitMap1;
  HitMap1 = new TH1D("HitMap1",title,1536,-0.5,1535.6);
  var+=">>HitMap1";
  myTree->Draw(var,myCuts);
  
  int max_hit=(int) HitMap1->GetMaximum();
  std::cout<<max_hit<<std::endl;
  //  gDirectory->Delete("PrjHitMap");
  TH1D *PrjHitMap;
  PrjHitMap = new TH1D("PrjHitMap","Projected Hit Map",100,0,max_hit);
  
  for(int strip=1;strip<=1535;strip++)
    {
      PrjHitMap->Fill(HitMap1->GetBinContent(strip));      
    }
  
  PrjHitMap->SetFillColor(5);
  //  int Integral = (int) PrjHitMap->GetSum();
  //  if (Integral!=0)  PrjHitMap->Scale(1./Integral);

  PrjHitMap->Draw();
  
  PrjHitMap->Fit("gaus","qsame");

  if(max_hit>=10)
    {
      TF1 *Gaussian = PrjHitMap->GetFunction("gaus");
      
      double sigma = Gaussian->GetParameter(2);
      double mean  = Gaussian->GetParameter(1);
      
      std::cout<<" Mean : "<<mean<<", Sigma "<<sigma<<std::endl;
      if(sigma>1) 
	{ 
	  for(int strip=1;strip<=1535;strip++)
	    {
	      int numHits = (int) HitMap1->GetBinContent(strip);
	      if(numHits==0)
		std::cout<<" WARNING!! Strip Number "<<strip-1<<" has 0 counts ****************************"<<std::endl;
	      
	      if(numHits > mean + 10.0*sigma)
		std::cout<<" WARNING!! Strip Number "<<strip-1<<" is above 10 sigma from the mean!! Counts : "<<numHits<<std::endl;
	    }
	}
    }

  myCanvas->Close();

  return PrjHitMap;
  
}

void PlotProjectHitMap(TString myCuts)
{
  TH1D *PHM = ProjectHitMap(myCuts);

  gDirectory->Delete("prjHitCan");
  TCanvas *prjHitCan;
  prjHitCan = new TCanvas("prjHitCan","Projected Hit Map",500,500);
  PHM->Draw();
}


TH1D *NumHitsLayer(TString myCuts)
{
  TString var="TkrNumHits";
  var.Prepend(TID);
  TString title=var;
  if(myCuts!="")
    {
      title+=" {";
      title+=myCuts;
      title+="}";
    }
  //  gDirectory->Delete("TkrNumHits");
  var+=">>TkrNumHits";
  TH1D *TkrNumHits = new TH1D("TkrNumHits",title,129,-0.5,128.5);

  myTree->Draw(var,myCuts);
  //  TH1D *TkrNumHits = (TH1D*) gDirectory->Get("TkrNumHits");
  TkrNumHits->SetFillColor(3);
  std::cout<<" Number of Selected Events: "<<TkrNumHits->GetEntries()
	   <<"("<<TkrNumHits->GetEntries()/NumberOfEvents*100.0<<"%)"<<std::endl;
  return   TkrNumHits;
}

void PlotNumHitsLayer(TString myCuts)
{
  gDirectory->Delete("NumHitsCan");
  TCanvas *NumHitsCan;
  NumHitsCan = new TCanvas("NumHitsCan","Num Hits",500,500);
  NumHitsCan->SetLogy();

  NumHitsLayer(myCuts);
}

long FindEvents(TString myCuts)
{
  return myTree->Scan("EventId",myCuts);
}

long FindBigEvents(int cut)
{
  int nx=11;  
  int ny=0;  
  TString myCut="";    
  
  for(int i=0; i<nx; i++)
    {
      TString anXlayer="X";
      anXlayer+=i;
      AddLayer(anXlayer);
      TString layerstring="Layer";
      layerstring +=anXlayer;
      myCut+=layerstring;
      myCut+=".TkrNumHits";
      myCut+=" >= ";
      myCut+= cut;
      if (i!=nx+ny-1) myCut += " || ";
    }


  for(int i=0; i<ny; i++)
    {
      TString anYlayer="Y";
      anYlayer+=i;
      AddLayer(anYlayer);
      TString layerstring="Layer";
      layerstring +=anYlayer;
      myCut+=layerstring;
      myCut+=".TkrNumHits";
      myCut+=" >= ";
      myCut+= cut;
      if (i!=ny-1) myCut += " || ";
    }
  
  return myTree->Scan("EventId",myCut);
}

TGraph *Correlation(TString var, TString VL1, TString VL2, TString myCuts)
{ 
  AddLayer(VL1);
  AddLayer(VL2);
  TString var0="Layer";
  TString var1="Layer";
  
  var0+=VL1;
  var1+=VL2;

  var0+=".";
  var1+=".";
  var0+=var;
  var1+=var;
  var0+=":";
  var0+=var1;
  
  myTree->Draw(var0,myCuts);
  TGraph *graph = (TGraph*) gPad->GetPrimitive("Graph");
  graph->SetMarkerStyle(25);
  graph->SetMarkerSize(0.2);
  graph->SetMarkerColor(2);
  graph->Draw("P");
  return graph;
}


TGraph *TimeHistory(TString VL1)
{

  AddLayer(VL1);
  TString var0="Layer";

  var0+=VL1;
  var0+=".TkrNumHits";
  var0+=":Header.EventId";
  myTree->Draw(var0);
  TGraph *graph= (TGraph*) gPad->GetPrimitive("Graph");
  graph->SetMarkerColor(2);
  graph->Draw("P");
  return graph;
  }

TGraph *Graph(TString VL1,  TString varA, TString varB, TString myCuts)
{ 
  AddLayer(VL1);
  TString var0="Layer";
  TString var1="Layer";
  
  var0+=VL1;
  var1+=VL1;

  var0+=".";
  var1+=".";
  var0+=varA;
  var1+=varB;
  var0+=":";
  var0+=var1;
  
  myTree->Draw(var0,myCuts);
  TGraph *graph = (TGraph*) gPad->GetPrimitive("Graph");
  graph->SetMarkerStyle(4);
  graph->SetMarkerSize(0.3);
  graph->SetMarkerColor(4);
  graph->Draw("P");
  return graph;
}


void PlotAllDigi(TString LV,TString myCuts)
{
  AddLayer(LV);
  PlotNumHitsLayer(myCuts);
  PlotHitMap(myCuts);
  PlotProjectHitMap(myCuts);
  PlotToT(myCuts);
  PlotGraph(LV,"ToT0","ToT1",myCuts);
}

void PlotCorrelation(TString var, TString VL1, TString VL2, TString myCuts)
{ 
  gDirectory->Delete("CorrCan");
  TCanvas *CorrCan;
  CorrCan = new TCanvas("CorrCan","Correlation between two Layers",500,500);
  Correlation(var,VL1,VL2,myCuts);
}

void PlotGraph(TString VL1,  TString varA, TString varB, TString myCuts)
{ 
  gDirectory->Delete("GraphCan");
  TCanvas *GraphCan;
  GraphCan = new TCanvas("GraphCan","Correlation between two variables",500,500);
  Graph(VL1,varA,varB,myCuts);
}

void PlotHitsVsTime(TString VL1)
{
  gDirectory->Delete("HistoryCan");
  TCanvas *HistoryCan;
  HistoryCan = new TCanvas("HistoryCan", "Hits vs Time",500,500);
  TimeHistory(VL1);
}

TCanvas *Report(TString LV,TString myCuts)
{
  //////////////////////////////////////////////////
  //  gDirectory->Delete("ReportCanvas");
  AddLayer(LV);
  TString Name="Off_line_Report ";
  TString Title="Off line Report ";
  Name+=myCuts;
  Title+=myCuts;
  TCanvas *ReportCanvas = new TCanvas("ReportCanvas","ReportCanvas",700,950);
  ReportCanvas->SetName(Name);
  ReportCanvas->SetTitle(Title);
  ReportCanvas->Divide(1,4);
  ReportCanvas->cd(1);
  TPad *top = (TPad*) gPad;

  ReportCanvas->cd(2);
  TPad *middleUp = (TPad*) gPad;

  ReportCanvas->cd(3);
  TPad *middleDown = (TPad*) gPad;

  ReportCanvas->cd(4);
  TPad *bottom = (TPad*) gPad;
  
  //  middleUp->Divide(2,1);
  middleDown->Divide(2,1);
  bottom->Divide(2,1);
  //////////////////////////////////////////////////

  TDatime date;
  int year    = date.GetYear();
  int month   = date.GetMonth();
  int day     = date.GetDay();
  int hour    = date.GetHour();
  int minute  = date.GetMinute();
  
  TString text;  
  TPaveText *pt = new TPaveText(0,0,1,1);
  pt->AddText(" Offline analysis report");
  TString text1=" Date : ";
  text1+=year;
  if(month<10)
    text1+=" / 0";
  else     
    text1+=" / 0";
  text1+=month;
  if(day<10)
    text1+=" / 0";
  else
    text1+=" / ";

  text1+=day;
  text1+="   -   ";
  text1+=hour;
  if (minute<10)
    text1+=":0";
  else 
    text1+=":";
  text1+=minute;
  pt->AddText(text1);


  text="Data File Name :";
  text+= myFile->GetName();  
  pt->AddText(text);

  pt->AddText(LayerName);

  /*
  text=" Layer = ";
  text+=LayerName;
  pt->AddText(text);
  */  


  pt->SetTextSize(.08);
  pt->SetTextColor(2);
  
  top->cd();
  pt->Draw();

 
  //  TH1D *SelectedHist;
  
  //  pt->AddText(myCuts);
  //////////////////////////////////////////////////
  //  long NumberOfEvents=myTree->GetEntries();

  TString text2=" TotalNumber of events :";
  text2+=NumberOfEvents;
  pt->AddText(text2);
  
  middleUp->cd();
  HitMap(myCuts);
  
  TH1D *HM = HitMap(myCuts);

  int nBins = HM->GetNbinsX();
  const TString sDeadBase = "dead strips:";
  TString sDead = sDeadBase;
  int nDead = 0;
  for ( int i=0; i<nBins; ++i ) {
      //      int value = (int)HM->GetBinContent(i+1);
      if ( (int)HM->GetBinContent(i+1) == 0 ) {
          sDead += " ";
          sDead += i;
          ++nDead;
      }
  }
  if ( nDead == 0 )
      sDead += " none";
  if ( sDead.Length() > 90 )
    {
      sDead = sDeadBase;
      sDead +=" ";
      sDead += nDead;
      sDead += " strips";
    }
  pt->AddText(sDead);

  //HM->SetName("HM");
  
  middleDown->cd(2);
  gPad->SetLogy();

  TH1D *NHL = NumHitsLayer(myCuts);
  
  long Sel_ToT1_ToT2_0_0 =  (long) NHL->GetEntries();

  TString text3 = " Number of selected events ( ";
  text3+=myCuts;
  text3+=" ) = ";
  text3+=Sel_ToT1_ToT2_0_0;
  double frac=1.0*Sel_ToT1_ToT2_0_0/(1.0*NumberOfEvents);
  char fractext[10];
  sprintf(fractext," (%.2g \%) ",frac*100.);
  text3+=fractext;

  pt->AddText(text3);
  std::cout<<text3<<std::endl;  
  
  middleDown->cd(1);
  //TH1D *PHM  = ProjectHitMap(myCuts);
  //middleDown->cd(1);
  //PHM->Draw();
  TGraph *TIMEGraph = TimeHistory(LV);
  
  
  bottom->cd(1);  
  gPad->SetLogy();
  TH1D *TOT0 = ToT_0(myCuts);
  
  
  bottom->cd(2);
  gPad->SetLogy();
  TH1D *TOT1  = ToT_1(myCuts);
  
  //ReportCanvas->cd();
  return ReportCanvas;
  //////////////////////////////////////////////////
}

void ReportAll(TString myCuts, TString opt) { 
  TString axis[2] = { "X", "Y" };
  for ( int i=0; i<18; ++i ) {
    for ( int iaxis=0; iaxis<2; ++iaxis ) {
      TString LV = axis[iaxis];
      LV += i;
      TString fileName = LV;
      fileName += ".";
      fileName += opt;
      Report(LV, myCuts)->Print(fileName);
    }
  }
}

void DumpReport(TString LV,char *filename)
{ 
  AddLayer(LV);
  TString PSFileName = "Report_";
  PSFileName+=LayerName;
  TString PSFileOpen=PSFileName;
  TString PSFileClose=PSFileName;
  PSFileOpen+=".ps(";
  PSFileClose+=".ps)";
  PSFileName+=".ps";

  Report(LV)->Print(PSFileOpen);
  
  Report(LV,"ToT0==0 && ToT1==0")->Print(PSFileName);
  Report(LV,"ToT0>0")->Print(PSFileName);
  Report(LV,"ToT1>0")->Print(PSFileName);
  Report(LV,"ToT0>0 && ToT1>0")->Print(PSFileName);
  Report(LV,"ToT0>0 || ToT1>0")->Print(PSFileName);
  Report(LV,"ToT0>249 && ToT1>249")->Print(PSFileName);
  Report(LV,"TkrNumHits >= 30")->Print(PSFileName);
  //  Report(LV,"TkrNumHits >= 60")->Print(PSFileName);
  //   Report(LV,"TkrNumHits >= 120")->Print(PSFileName);
  
  TCanvas *CorrCan = new TCanvas("CorrCan","CorrelationsBetween Layer",700,700);
  CorrCan->Divide(2,3);
  CorrCan->cd(1);
  Correlation("TkrNumHits","X0","X1","LayerX0.TkrNumHits>0");
  CorrCan->cd(2);
  Correlation("TkrNumHits","X0","X2","LayerX0.TkrNumHits>0");
  CorrCan->cd(3);
  Correlation("TkrNumHits","X0","X3","LayerX0.TkrNumHits>0");
  CorrCan->cd(4);
  Correlation("TkrNumHits","X1","X2","LayerX1.TkrNumHits>0");
  CorrCan->cd(5);
  Correlation("TkrNumHits","X1","X3","LayerX1.TkrNumHits>0");
  CorrCan->cd(6);
  Correlation("TkrNumHits","X2","X3","LayerX2.TkrNumHits>0");
  /*
  Correlation("TkrNumHits","X0","X1","LayerX0.TkrNumHits>32");
  CorrCan->cd(2);
  Correlation("TkrNumHits","X0","X2","LayerX0.TkrNumHits>32");
  CorrCan->cd(3);
  Correlation("TkrNumHits","X0","X3","LayerX0.TkrNumHits>32");
  CorrCan->cd(4);
  Correlation("TkrNumHits","X1","X2","LayerX1.TkrNumHits>32");
  CorrCan->cd(5);
  Correlation("TkrNumHits","X1","X3","LayerX1.TkrNumHits>32");
  CorrCan->cd(6);
  Correlation("TkrNumHits","X2","X3","LayerX2.TkrNumHits>32");
  */
  CorrCan->Print(PSFileClose);
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
  gDirectory->Delete("CanVtx");
  CanVtx = new TCanvas("CanVtx","CanVtx",1000,300);
  CanVtx->Divide(3,1);
  CanVtx->cd(1);
  myTree->Draw("TkrVtxX",myCuts);
  CanVtx->cd(2);
  myTree->Draw("TkrVtxY",myCuts);
  CanVtx->cd(3);
  myTree->Draw("TkrVtxZ",myCuts);
}

TH2D *ClusHit(TString myCuts)
{
  TString var0="TkrNumClus";
  TString var1="TkrTotalNumHits";
  TString var=var0;
  var+=":";
  var+=var1;
  TString title=var;
  if(myCuts!="")
    {
      title+=" {";
      title+=myCuts;
      title+="}";
    }
//  gDirectory->Delete("NumClus_NumHits");
  TH2D *NumClus_NumHits;
  NumClus_NumHits = new TH2D("NumClus_NumHits",var,30,0,15,30,0,15);
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
  return NumClus_NumHits;
}

void PlotClusHit(TString myCuts)
{
  gDirectory->Delete("ClusCan");
  TCanvas *ClusCan;
  ClusCan = new TCanvas("ClusCan","Clusters vs Hits",700,700);
  ClusHit(myCuts);
}

void PlotRecon(TString myCuts)
{
  gDirectory->Delete("RecCan");
  TCanvas *RecCan;
  RecCan = new TCanvas("RecCan","Recontruction Canvas",700,700);
  RecCan->Divide(2,2);
  RecCan->cd(1);
  myTree->Draw("Theta",myCuts);
  RecCan->cd(2);
  myTree->Draw("Phi",myCuts);
  RecCan->cd(3);
  myTree->Draw("ThetaXZ",myCuts);
  RecCan->cd(4);
  myTree->Draw("ThetaYZ",myCuts);
  RecCan->cd(1);
}


void CoincidenceBigEvents(int NH)
{
  int NumXlayers = 11;
  int NumYlayers = 0;  
  int NumLayers = NumXlayers + NumYlayers;
  std::cout<<NumLayers<<std::endl;
  int *NumHits = new int[NumLayers];
  
  long EventId;
  myTree->SetBranchAddress("EventId", &EventId);
  for(int nx = 0; nx < NumXlayers; nx++)
    {
      TString anXLayer = "X";
      anXLayer += nx;
      AddLayer(anXLayer);
      TString myVar = "Layer";
      myVar+=anXLayer;
      myVar+=".TkrNumHits";
      myTree->SetBranchAddress(myVar,&NumHits[nx]);
    }
  
  for(int ny = 0; ny < NumYlayers; ny++)
    {  
      TString anYLayer = "Y";
      anYLayer += ny;
      AddLayer(anYLayer);
      TString myVar = "Layer";
      myVar+=anYLayer;
      myVar+=".TkrNumHits";
      myTree->SetBranchAddress(myVar,&NumHits[NumXlayers+ny]);
    }

  int N = myTree->GetEntries();
  
  TString HistTitle="Fraction of events vs Coincidence layers (NumHits > ";
  HistTitle += NH;
  HistTitle +=")";
  
  TH1D *myHist = new TH1D("Events",HistTitle,NumLayers,.5,NumLayers+.5);

  for(int i=0; i<N; i++)
    {
      myTree->GetEntry(i);
      
      int myLayers = 0;
      for(int l = 0; l<NumLayers; l++)
	{
	  if(NumHits[l] > NH) myLayers++;
	}
      if(myLayers>0) 
	{
	  myHist->Fill(myLayers);
	  std::cout<<EventId<<std::endl;
	}
    }
  myHist->Scale(1./N);
  myHist->Draw();
  
}

/*
void EventDisplay(long EventId)
{
  int NumXlayers = 11;
  int NumYlayers = 8;  

  int NumLayers = NumXlayers + NumYlayers;
  
  std::cout<<NumLayers<<std::endl;
  
  TString myCuts = "EventId==";
  myCuts+=EventId;
  
  float XStripId[128*18];
  float XZposition[128*18];

  float YStripId[128*18];
  float YZposition[128*18];
  
  float XLayersPositions[18] = {-16.8, -13.2, -23.3, -19.7, -29.8, -26.2, -36.2, -32.6, -42.7, -39.1, -55.7, -45.6, -55.7, -52.1, -62.2, -58.6, -68.7, -65.1};
  float YLayersPositions[18] = {16.8, 13.2, 23.3, 19.7, 29.8, 26.2, 36.2, 32.6, 42.7, 39.1, 55.7, 45.6, 55.7, 52.1, 62.2, 58.6, 68.7, 65.1};
  
  
  //  std::cout<<TotalNumHits<<std::endl;
  int idx=0;

  // X Layers:

  for(int nx = 0; nx < NumXlayers; nx++)
    {
      TString anXLayer = "X";
      anXLayer += nx;
      AddLayer(anXLayer);

      gDirectory->Delete("HitMap");
      TH1D *HitMapHisto = HitMap(myCuts);
      int NumHits = HitMapHisto->GetEntries();
      if(NumHits>0)
	{
	  int strip=1;
	  int hitperlayer=0;
	  while(strip <= 1536 && hitperlayer <= NumHits)
	    {
	      if(HitMapHisto->GetBinContent(strip)>0)
		{
		  XStripId[idx]  = (float) HitMapHisto->GetBinCenter(strip);
		  XZposition[idx]= (float) XLayersPositions[nx];
		  idx++;
		  hitperlayer++;
		}
	      strip++;
	    }
	}
      
      
    }
  TGraph *gX = new TGraph(idx,XStripId,XZposition);
  gX->SetTitle("X Layers");
  
  // Y Layers:
  int idy=0;
  

  for(int ny = 0; ny < NumYlayers; ny++)
    {
      TString anYLayer = "Y";
      anYLayer += ny;
      AddLayer(anYLayer);

      gDirectory->Delete("HitMap");
      TH1D *HitMapHisto = HitMap(myCuts);
      
      int NumHits = HitMapHisto->GetEntries();
      
      if(NumHits>0)
	{
	  int strip=1;
	  int hitperlayer=0;
	  while(strip <= 1536 && hitperlayer <= NumHits)
	    {
	      if(HitMapHisto->GetBinContent(strip)>0)
		{
		  YStripId[idy]  = (float) HitMapHisto->GetBinCenter(strip);
		  YZposition[idy]= (float) YLayersPositions[ny];
		  idy++;
		  hitperlayer++;
		}
	      strip++;
	    }
	}
      
      
    }
  
  TGraph *gY = new TGraph(idy,YStripId,YZposition);
  gY->SetTitle("Y Layers");

  //////////////////////////////////////////////////
  //                     Dispaly
  TCanvas *EventDisplayCanvas = new TCanvas("EventDisplayCanvas","EventDisplayCanvas",800,800);
  EventDisplayCanvas->Divide(1,2);
  EventDisplayCanvas->cd(2);
  gX->Draw("ap*");
  EventDisplayCanvas->cd(1);
  gY->Draw("ap*");
}
*/

void AllLayers(int cut=0)
{
  int nx=11;  
  int ny=0;  
  TString myCut="";    

  for(int i=0; i<nx; i++)
    {
      TString anXlayer="X";
      anXlayer+=i;
      AddLayer(anXlayer);
      TString layerstring="Layer";
      layerstring +=anXlayer;
      myCut+=layerstring;
      myCut+=".TkrNumHits";
      myCut+=" >= ";
      myCut+= cut;
      if (i!=nx+ny-1) myCut += " || ";
    }


  for(int i=0; i<ny; i++)
    {
      TString anYlayer="Y";
      anYlayer+=i;
      AddLayer(anYlayer);
      TString layerstring="Layer";
      layerstring +=anYlayer;
      myCut+=layerstring;
      myCut+=".TkrNumHits";
      myCut+=" >= ";
      myCut+= cut;
      if (i!=ny-1) myCut += " || ";
    }

  std::cout<<myCut<<std::endl;

  TH1D* h0;

  for(int i=0; i < nx+ny; i++)
    {
      if(i<nx)
	{
	  TString anXlayer="X";
	  anXlayer+=i;
	  AddLayer(anXlayer);
	  std::cout<<anXlayer<<std::endl;
	}
      else //if(i<ny)
	{
	  TString anYlayer="Y";
	  int ii = (i-nx);
	  anYlayer+=ii;
	  AddLayer(anYlayer);
	  std::cout<<anYlayer<<std::endl;
	}
      
      if(i==0) h0 = NumHitsLayer(myCut);    
      if(i>0) h0->Add(NumHitsLayer(myCut)); 
    }
  h0->Draw();
}

void PlotTotalTkrNumHits()
{
  TH1D *TotalNumHits = new TH1D("TotalNumHits","TotalNumHits",128,0,128);
  myTree->Draw("TkrTotalNumHits>>TotalNumHits");
  TotalNumHits->Draw();
}

void AllLayers01(int cut)
{
  int nx=11;  
  int ny=0;  
  TString myCut="";    

  for(int i=0; i<nx; i++)
    {
      TString anXlayer="X";
      anXlayer+=i;
      AddLayer(anXlayer);
      TString layerstring="Layer";
      layerstring +=anXlayer;
      myCut+=layerstring;
      myCut+=".TkrNumHits";
      myCut+=" <= ";
      myCut+= cut;
      if (i!=nx+ny-1) myCut += " && ";
    }


  for(int i=0; i<ny; i++)
    {
      TString anYlayer="Y";
      anYlayer+=i;
      AddLayer(anYlayer);
      TString layerstring="Layer";
      layerstring +=anYlayer;
      myCut+=layerstring;
      myCut+=".TkrNumHits";
      myCut+=" <=  ";
      myCut+= cut;
      if (i!=ny-1) myCut += " && ";
    }

  std::cout<<myCut<<std::endl;

  TH1D* h0;

  for(int i=0; i < nx+ny; i++)
    {
      if(i<nx)
	{
	  TString anXlayer="X";
	  anXlayer+=i;
	  AddLayer(anXlayer);
	  std::cout<<anXlayer<<std::endl;
	}
      else //if(i<ny)
	{
	  TString anYlayer="Y";
	  int ii = (i-nx);
	  anYlayer+=ii;
	  AddLayer(anYlayer);
	  std::cout<<anYlayer<<std::endl;
	}
      
      if(i==0) h0 = NumHitsLayer(myCut);    
      if(i>0) h0->Add(NumHitsLayer(myCut)); 
    }
  h0->Draw();
}
