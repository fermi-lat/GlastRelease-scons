#include <iomanip>
#include <iostream>
#include <fstream>
#include <string>

#include "TGraph.h"
#include "TTree.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TString.h"
#include "TText.h"
#include "TMap.h"
#include "TObjString.h"
#include "TLine.h"
#include "TF1.h"
#include "TGraph.h"

TFile* myFile;
TMap* myGeometry;

void EventDisplay(TString="MyRootFile.root", int=0, int=-1);
bool isGeometry();
void loadGeometry(TString filename="src/LeaningTower/geometry/stack2geometry.txt");
double getGeometry(TString layer);

bool isGeometry() { return myGeometry && !myGeometry->IsEmpty(); }

void loadGeometry(TString filename) {
    if ( isGeometry() )
        myGeometry->DeleteAll();
    if ( !myGeometry ) {
        std::cout << "making NEW geometry map" << std::endl;
        myGeometry = new TMap();
    }
    std::ifstream fin(filename);
    if ( !fin.is_open() )
        std::cerr << "ERROR:  File " << filename << " couldn't be opened!" << std::endl;
    TString layer;
    TString value;
    while ( !fin.eof() ) {
        fin >> layer >> value;
        if ( fin.good() )
            myGeometry->Add(new TObjString(layer), new TObjString(value));
    }
    fin.close();
    myGeometry->Print();
}

double getGeometry(TString layer) {
    if ( !isGeometry() ) {
        std::cerr << "ERROR:  Geometry is empty!  Maybe you forgot to load it?" << std::endl;
        return 1E15;
    }
    if ( !myGeometry->FindObject(layer) ) {
        std::cerr << "ERROR: Layer " << layer << " not found in geometry map!" << std::endl;
        return 1E15;
    }
    return atof(((TObjString*)myGeometry->GetValue(layer))->GetName());
}

//////////////////////////////////////////////////
double returnCoordinate(int StripNumber){
  bool DEBUG = false;
  // Check that the strip number is within the allowed range.
  if ((StripNumber < 0) || (StripNumber > 1535)){
    std::cout << "WARNING: the strip number must be included in the range [0:1535]." << std::endl;
    std::cout << "returnCoordinate() returning -1.0" << std::endl;
    return -1.0;
  }
  // Define some constants (to be moved out of the function?)
  // All dimension in cm.
  double EDGE_WIDTH        = 0.1000;
  double STRIP_PITCH       = 0.0228;
  double LADDER_SEPARATION = 0.0200;

  double Coordinate = EDGE_WIDTH + STRIP_PITCH*StripNumber +
  (LADDER_SEPARATION + 2*EDGE_WIDTH - STRIP_PITCH)*(int)(StripNumber/384);
  if (DEBUG)
    {
      std::cout << "returnCoordinate() returning " << Coordinate <<
	"(strip number = " << StripNumber << ")" << std::endl;
    }
  return Coordinate;
}

bool checkActiveArea(double ParallelCoordinate, double NormalCoordinate, double BorderWidth){
  bool DEBUG = false;
  // Define some constants (to be moved out of the function?)
  // All dimension in cm.
  double EDGE_WIDTH        = 0.1000;
  double WAFER_WIDTH       = 8.9500;
  double LADDER_SEPARATION = 0.0200;
  // Loop over the ladders.
  for (int i=0; i<4; i++){
    if ((NormalCoordinate > (EDGE_WIDTH + (WAFER_WIDTH + LADDER_SEPARATION)*i + BorderWidth)) &&
	(NormalCoordinate < (WAFER_WIDTH*(i+1) + LADDER_SEPARATION*i - BorderWidth))){
      for (int j=0; j<4; j++){
	if ((ParallelCoordinate > (EDGE_WIDTH + WAFER_WIDTH*i + BorderWidth)) &&
	    (ParallelCoordinate < (WAFER_WIDTH*(i+1) - BorderWidth))){
	  if (DEBUG){
	    std::cout << "checkActiveArea(): ladder " << i << ", wafer " << j << " hit." << std::endl;
	  }
	  return 1;
	}
      }
    }
  }
  if (DEBUG){
    std::cout << "checkActiveArea(): no hit on any of the wafers." << std::endl;
  }
  return 0;
}

//////////////////////////////////////////////////
void EventDisplay(TString filename, int firstEvent, int lastEvent) 
{
  
  bool DISPLAY=true;
  int gap=5;

  if ( !myFile ) {
      std::cout << "opening new myFile" << std::endl;
      myFile = new TFile(filename);
  }
  TTree *myTree = (TTree*)gDirectory->Get("Header");
  int numEvents = myTree->GetEntries();
  std::cout << "Number of Events: " << numEvents << std::endl;

  TCanvas *myCanvas = new TCanvas("c1","Event Display", 700, 800);

  Int_t EventId;
  Int_t RunId;
  Int_t TkrTotalNumHits;  
  myTree->SetBranchAddress("EventId",&EventId);
  myTree->SetBranchAddress("RunId",&RunId);
  myTree->SetBranchAddress("TkrTotalNumHits",&TkrTotalNumHits);
  //int iEvent=firstEvent;
  int IneffX[18];
  TString LayerLabelX[18];
  TString LayerLabelY[18];
  int IneffY[18];
  int EventsX[18];
  int EventsY[18];

  for(int i = 0; i<18; i++)
    {
      EventsX[i]= 0;
      EventsY[i]= 0;
      IneffX[i] = 0;
      IneffY[i] = 0;
    }
  
  int NumberOfLayersX,NumberOfLayersY;
  
  for ( int iEvent=firstEvent; iEvent<numEvents; ++iEvent ) 
    {
      std::cout<<iEvent<<std::endl;
      myTree->GetEntry(iEvent);
      
      std::cout << "EventId: " << EventId
		<< " RunId: " << RunId
		<< " TkrTotalNumHits: "  << TkrTotalNumHits
		<< std::endl;

      if ( !isGeometry() ) {
	std::cerr << "ERROR:  Geometry is empty! Maybe you forgot to load it?" << std::endl;
	return;
      }
      TMapIter ti(myGeometry);
      TObjString* key;
      TTree* layer;
      double x[4] = { 0, 40, 40,0 };
      double y[4] = { -65, -65, 65, 65 };
      
      TGraph *tg;
      TText *text1,*text2;
      TLine *layerLine;
      
      if(DISPLAY)
	{
	  tg = new TGraph(4, x, y);
	  tg->GetXaxis()->SetTitle("strip Id");
	  tg->GetYaxis()->SetTitle("position in stack/cm");
	  tg->Draw("AP");
	  TString title = "EventId: ";
	  title+= EventId;
	  title+=" TkrTotalNumHits: ";
	  title+=TkrTotalNumHits;
	  text1 = new TText(0, 70, title);
	  text1->SetTextAlign();
	  text1->SetTextSize(0.04);
	  text1->Draw();
	}
      
      //////////////////////////////////////////////////
      int recon=1;
      double *Xx = new double[TkrTotalNumHits];
      double *Zx = new double[TkrTotalNumHits];
      
      double *Xy = new double[TkrTotalNumHits];
      double *Zy = new double[TkrTotalNumHits];
      
      int hitsX=0;
      int hitsY=0;
      
      double ThetaX,CosThX,SinThX;
      double StackPositionsX[18];
      double StackPositionsY[18];
      int    NumHitsLayerX[18];
      int    NumHitsLayerY[18];
      NumberOfLayersX=0;
      NumberOfLayersY=0;
      
      
      //////////////////////////////////////////////////
      
      while ( key = (TObjString*)ti.Next() ) 
	{
	  TString VL = key->String();
	  double pos = getGeometry(VL);
	  TString LayerName = "Layer";
	  LayerName+=VL;
	  layer = (TTree*)myFile->Get(LayerName);
	  Int_t TkrNumHits;
	  Int_t TkrHits[128];
	  layer->SetBranchAddress("TkrNumHits",&TkrNumHits);
	  layer->SetBranchAddress("TkrHits",TkrHits);
	  layer->GetEntry(iEvent);
	  
	  if(DISPLAY)
	    { 
	      layerLine = new TLine(returnCoordinate(0), pos, returnCoordinate(1535),pos);
	      layerLine->Draw();
	      TString title = VL;
	      title+=" (";
	      title+=TkrNumHits;
	      title+=")";
	      text2 = new TText(40, pos, title);
	      text2->SetTextAlign(02);
	      text2->SetTextSize(0.02);
	      text2->Draw();
	    }
	  

	  //std::cout<<"TkrNumHits : "<<VL<<" "<<TkrNumHits<<std::endl;
	  
	  //////////////////////////////////////////////////
	  if(pos>0)
	    {
	      NumHitsLayerY[NumberOfLayersY]     = TkrNumHits;
	      LayerLabelY[NumberOfLayersY] = VL;
	      StackPositionsY[NumberOfLayersY++] = pos;		
	      
	    }
	  else
	    {
	      LayerLabelX[NumberOfLayersX] = VL;
	      NumHitsLayerX[NumberOfLayersX]     = TkrNumHits;
	      StackPositionsX[NumberOfLayersX++] = pos;
	    }	    
	  
	  
	  if ( TkrTotalNumHits == 0) recon=0;
	  else if ( TkrNumHits > 0) 
	    {
	      Double_t xpos[128];
	      Double_t zpos[128];
	      
	      for ( int i=0; i<TkrNumHits; ++i ) 
		{
		  xpos[i] = returnCoordinate(TkrHits[i]);
		  zpos[i] = pos;
		  
		  if(pos>0) 
		    {
		      Zy[hitsY]=zpos[i];
		      if(i==0) 
			{
			  Xy[hitsY] = xpos[i];
			}
		      else if(TkrHits[i]<= TkrHits[i-1]+gap && TkrHits[i] >= TkrHits[i-1]-gap) 
			{
			  Xy[hitsY] = (xpos[i]+xpos[i-1])/2.;
			}		      
		      else recon=0;
		    }
		  else
		    {
		      Zx[hitsX]=zpos[i];
		      if(i==0) 
			{
			  Xx[hitsX] = xpos[i];
			}
		      else if(TkrHits[i]<= TkrHits[i-1]+gap && TkrHits[i]>= TkrHits[i-1]-gap) 
			{
			  Xx[hitsX] = (xpos[i]+xpos[i-1])/2.;
			}
		      else recon=0;
		    }
		}
	      
	      //////////////////////////////////////////////////
	      if(pos>0) 
		{
		  hitsY++;
		}
	      else
		{
		  hitsX++;
		}
	      //////////////////////////////////////////////////
	      
	      if(DISPLAY) 
		{ 
		  TGraph *hitsGraph = new TGraph(TkrNumHits, xpos, zpos);
		  hitsGraph->Draw("*");
		}
	    }
	}
      
      if ( hitsY*hitsX == 0) recon=0;
      //      std::cout<<"hits X = "<<hitsX<<" hits Y = "<<hitsY<<" Recon = "<<recon<<std::endl;

      TGraph *Xtrack,*Ytrack;
      TF1 *fX,*fY; 

      if(recon)
	{
	  double Xextrapolated,Yextrapolated;
	  
	  Xtrack = new TGraph(hitsX,Xx,Zx);
	  fX = new TF1("fX","pol1",0,40);
	  
	  Ytrack = new TGraph(hitsY,Xy,Zy);
	  fY = new TF1("fY","pol1",0,40);
	  
	  
	  if(DISPLAY)
	    {
	      Xtrack->SetMarkerStyle(24);
	      Xtrack->SetMarkerColor(4); 
	      //	    Xtrack->Draw();
	      fX->SetLineColor(2);
	      fX->SetLineStyle(2);
	      Xtrack->Draw("p");
	      Ytrack->SetMarkerStyle(24);
	      Ytrack->SetMarkerColor(4);
	      fY->SetLineColor(3);
	      fY->SetLineStyle(2);
	      Ytrack->Draw("p");
	    }
	  
	  Xtrack->Fit("fX","QR");	  
	  Ytrack->Fit("fY","QR");

	  ///////// Check for Inefficiences
	  double Ax = fX->GetParameter(0);
	  double Bx = fX->GetParameter(1);
	  
	  double Ay = fY->GetParameter(0);
	  double By = fY->GetParameter(1);
	  
	  for(int layerY = 0; layerY < NumberOfLayersY; layerY++)
	    {
	      double Zlayer = StackPositionsY[layerY];
	      Yextrapolated = (Zlayer - Ay)/By;
	      Xextrapolated = (Zlayer - Ax)/Bx;
	      
	      bool ActiveRegion=false;
	      if(checkActiveArea(Xextrapolated,Yextrapolated,1.0))
		{
		  ActiveRegion=true;
		  EventsY[layerY]++;
		}
	      if(ActiveRegion && NumHitsLayerY[layerY]==0) 
		{
		  IneffY[layerY]++;
		}
	      
	    }
	  
	  for(int layerX = 0; layerX < NumberOfLayersX;layerX++)
	    {
	      double Zlayer = StackPositionsX[layerX];
	      Yextrapolated = (Zlayer - Ay)/By;
	      Xextrapolated = (Zlayer - Ax)/Bx;
	      bool ActiveRegion=false;
	      if(checkActiveArea(Yextrapolated,Xextrapolated,1.0)) 
		{
		  ActiveRegion=true;
		  EventsX[layerX]++;
		}
	      if(ActiveRegion && NumHitsLayerX[layerX]==0) 
		{
		  IneffX[layerX]++;
		}
	      
	    }
	  
	  
	}
      
      
      
      
      //////////////////////////////////////////////////
      
      
      //////////////////////////////////////////////////
      
      if(DISPLAY)  myCanvas->Update();

      if(recon)
	{
	  for(int layerY = 0; layerY < NumberOfLayersY; layerY++)
	    if(EventsY[layerY]>0)
	      std::cout << "Inefficiences  "<<LayerLabelY[layerY]<<" = "<<IneffY[layerY]<<" "<<EventsY[layerY]<<" "<<(float) IneffY[layerY]/EventsY[layerY]<<std::endl;
	  
	  for(int layerX = 0; layerX < NumberOfLayersX;layerX++)
	    if(EventsX[layerX]>0)
	      std::cout << "Inefficiences  "<<LayerLabelX[layerX]<<" = "<<IneffX[layerX]<<" "<<EventsX[layerX]<<" "<<(float) IneffX[layerX]/EventsX[layerX]<<std::endl;
	}
      
      if(lastEvent<iEvent)  
	{    
	  std::cout << "Hit return to continue (q to quit): ";
	  char c = getchar();
	  //////////////////////////////////////////////////
	  if ( c == 'q' )
	    return;
	  else 
	    {
	      if(DISPLAY)
		{
		  delete tg;
		  delete text1;
		  delete text2;
		  delete layerLine;
		}
	      delete Xx;
	      delete Zx;
	      delete Xy;
	      delete Zy;
	      
	      if(recon)
		{
		  delete Xtrack;
		  delete Ytrack;
		  delete fX;
		  delete fY; 
		}
	      
	    }
	}
      if(iEvent == lastEvent)
	return;
    }
}

