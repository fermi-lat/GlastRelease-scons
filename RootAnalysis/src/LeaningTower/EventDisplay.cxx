#include "TGraph.h"
#include "TGraphErrors.h"
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
#include "TStyle.h"
#include "TSystem.h"

#include <vector>
#include "Tracker.h"
#include "Layer.h"
#include "Event.h"

Event *myEvent;
Tracker *tracker;
TCanvas *EventDisplayC;

TText LabelNumHits[36];

TGraph *XTrack;
TGraph *YTrack;

TGraph *XClusters;
TGraph *YClusters;

int EventId, RunId,TkrTotalNumHits;

void InitializeED(TString filename = "MyRootFile.root")
{
  //  Initialize(filename);
  gStyle->SetCanvasColor(10);
  tracker = new Tracker();

  tracker->loadGeometry(gSystem->ExpandPathName("$ROOTANALYSISROOT/src/LeaningTower/geometry/stack2geometry.txt"));
  
  myEvent = new Event(filename,(TMap *) tracker->GetGeometry());
  
  EventDisplayC = new TCanvas("EventDisplayC"  ,"EventDisplayC",1100,800);  
  tracker->Display();
  
  //////////////////////////////////////////////////

  XTrack = new TGraph();
  YTrack = new TGraph();
  XTrack->SetMarkerStyle(7); 
  YTrack->SetMarkerStyle(7);
  //////////////////////////////////////////////////
  XClusters = new TGraph();
  YClusters = new TGraph();

  XClusters->SetMarkerStyle(4); 
  YClusters->SetMarkerStyle(4);
  XClusters->SetMarkerColor(3); 
  YClusters->SetMarkerColor(3);
}

void DisplayEvent(int NumEvent=0)
{

  myEvent->Go(NumEvent);
  
  int TkrTotalNumHits = myEvent->GetTkrTotalNumHits();


  std::cout<<" TkrTotalNumHits = "<<TkrTotalNumHits<<std::endl;
  std::cout<<" EventId         = "<<myEvent->GetEventId()<<std::endl;
  std::cout<<" RunId           = "<<myEvent->GetRunId()<<std::endl;

  
  TMap *myGeometry = tracker->GetGeometry();
  TMapIter ti(myGeometry);
  TObjString* key;
  int i=0;
  
  
  double *XHitXlayers = new double[TkrTotalNumHits];
  double *ZHitXlayers = new double[TkrTotalNumHits];
  double *XHitYlayers = new double[TkrTotalNumHits];
  double *ZHitYlayers = new double[TkrTotalNumHits];

  int Xhits = 0;
  int Yhits = 0;

  std::vector<double> XClusterXlayer;
  std::vector<double> ZClusterXlayer;
  std::vector<double> XClusterYlayer;
  std::vector<double> ZClusterYlayer;
  
  while ( ( key = (TObjString*)ti.Next() ) ) 
    {
      Layer *aLayer = ((Layer*) myGeometry->GetValue(key));
      
      TString LayerName=key->String();

      double height = aLayer->GetHeight();
      
      int LayerNumHits = myEvent->GetLayerNumHits(LayerName);
      int *LayerHits   = myEvent->GetLayerHits(LayerName);
      //      std::cout<<"GetLayerNumHits ("<<LayerName<<") = "<<LayerNumHits<<" "<<height<<std::endl;
      
      TString title = "(";
      title        += LayerNumHits;
      title        += ")";

      LabelNumHits[i].SetText(40.0, height, title);
      LabelNumHits[i].SetTextSize(0.02);
      if(aLayer->IsX()) 
	EventDisplayC->cd(1);      
      else
	EventDisplayC->cd(2);
      LabelNumHits[i].Draw();
      i++;      

      if(LayerNumHits>0)
	{
	  //////////////////////////////////////////////////
	  //                   Clustering
	  std::vector<double> ClusterLayer = myEvent->GetClusters(LayerName);
	  int ClusterLayerSize = ClusterLayer.size();
	  //////////////////////////////////////////////////
	  //	  std::cout<<LayerName<<" hits = "<<LayerNumHits<<" clusters = "<<ClusterLayerSize <<std::endl;
	  
	  if(aLayer->IsX()) 
	    {
	      EventDisplayC->cd(1);
	      // Fill the Hits arrays...
	      
	      for(int j = 0; j < LayerNumHits; j++)
		{
		  XHitXlayers[Xhits] = aLayer->GetCoordinate(LayerHits[j]);
		  ZHitXlayers[Xhits] = height;
		  Xhits++;
		}
	      
	      // Fill the clusters
	      
	      for (int j=0; j<ClusterLayerSize;j++)
		{
		  XClusterXlayer.push_back(ClusterLayer[j]);
		  ZClusterXlayer.push_back(height);
		}
	      
	      //      std::cout<<"GetLayerNumHits ("<<LayerName<<") = "<<LayerNumHits<<" "<<height<<std::endl;
	    } 
	  else
	    { 
	      EventDisplayC->cd(2);
	      // Fill the Hits arrays...
	      for(int j=0; j<LayerNumHits; j++)
		{
		  XHitYlayers[Yhits] = aLayer->GetCoordinate(LayerHits[j]);
		  ZHitYlayers[Yhits] = height;
		  Yhits++;
		}
	      
	      // Fill the clusters
	      
	      for (int j=0; j<ClusterLayerSize;j++)
		{
		  XClusterYlayer.push_back(ClusterLayer[j]);
		  ZClusterYlayer.push_back(height);
		}
	      
	    }
	}
    }
  
  // Display the track(s):



  
  XTrack->Set(Xhits);
  YTrack->Set(Yhits);  
  
  for (int i =0;i<Xhits;i++)
    {
      XTrack->SetPoint(i,XHitXlayers[i],ZHitXlayers[i]);
      
    }
  for (int i =0;i<Yhits;i++)
    {
      YTrack->SetPoint(i,XHitYlayers[i],ZHitYlayers[i]);
    }


  EventDisplayC->cd(1);
  XTrack->Draw("P");
  EventDisplayC->cd(2);
  YTrack->Draw("P");
  
  // Display the Cluster(s)

  int NumClusX = XClusterXlayer.size();
  int NumClusY = XClusterYlayer.size();  
  
  XClusters->Set(NumClusX);
  YClusters->Set(NumClusY);

  for (int i =0;i<NumClusX;i++)
    {
      XClusters->SetPoint(i,XClusterXlayer[i],ZClusterXlayer[i]);
    }
  for (int i =0;i<NumClusY;i++)
    {
      YClusters->SetPoint(i,XClusterYlayer[i],ZClusterYlayer[i]);
    }

  EventDisplayC->cd(1);
  XClusters->Draw("P");
  gDirectory->Delete("fun");
  //  TF1 *fun = new TF1("fun","pol1",0,40);
  TF1 fun("fun","pol1",0,40);
  //  fun->SetParameters(0,-1000);
  //  fun->SetParLimits(0,-1000,1000);
  //  fun->SetParLimits(1,-1000,1000);
  
  //  fun->SetLineStyle(2);
  //  fun->SetLineWidth(0);
  //  fun->SetLineColor(2);
  fun.SetLineStyle(2);
  fun.SetLineWidth(0);
  fun.SetLineColor(2);

  if(NumClusX>1)
    {
      XClusters->Fit("fun","R");
    }
  else
    {
        TObject* dummy = XClusters->FindObject("fun");
        if ( dummy ) dummy->Delete();
    }
  
  EventDisplayC->cd(2);
  YClusters->Draw("P");
  if(NumClusY>1)
    {
      YClusters->Fit("fun","R");
    } 
  else
    {
        TObject* dummy = YClusters->FindObject("fun");
        if ( dummy ) dummy->Delete();
    }
  EventDisplayC->Update();
}


void IneffAnalysis(int LastEvent, TString LV="All")
{
  std::vector<double> XInActiveArea;
  std::vector<double> YInActiveArea;

  std::vector<double> XNotInActiveArea;
  std::vector<double> YNotInActiveArea;  

  std::vector<double> YMissedHits;
  std::vector<double> XMissedHits;

  gDirectory->Delete("XResiduals");
  gDirectory->Delete("YResiduals");

  TH1D *XResiduals = new TH1D("XResiduals","XResiduals",100,-1,1);
  TH1D *YResiduals = new TH1D("YResiduals","YResiduals",100,-1,1);

  TMap *myGeometry = myGeometry = tracker->GetGeometry();
  TMapIter ti(myGeometry);
  TObjString* key;

  for(int ev = 1;ev < LastEvent; ev++)
    { 
      
      myEvent->Go(ev);
      
      int TkrTotalNumHits = myEvent->GetTkrTotalNumHits();
      if(ev%1000==0)
	{
	  std::cout<<" TkrTotalNumHits = "<<TkrTotalNumHits<<std::endl;
	  std::cout<<" EventId         = "<<myEvent->GetEventId()<<std::endl;
	  std::cout<<" RunId           = "<<myEvent->GetRunId()<<std::endl;
	}
      
      ti.Reset();
            
      std::vector<double> XClusterXlayer;
      std::vector<double> ZClusterXlayer;
      std::vector<double> XClusterYlayer;
      std::vector<double> ZClusterYlayer;

      while ( ( key = (TObjString*)ti.Next() ) ) 
	{
	  Layer *aLayer = ((Layer*) myGeometry->GetValue(key));
       	  TString LayerName=key->String();
	  
	  double height    = aLayer->GetHeight();
	  int LayerNumHits = myEvent->GetLayerNumHits(LayerName);
	  
	  if(LayerNumHits>0)
	    {
	      //////////////////////////////////////////////////
	      //                   Clustering
	      std::vector<double> ClusterLayer = myEvent->GetClusters(LayerName);
	      int ClusterLayerSize             = ClusterLayer.size();
	      //////////////////////////////////////////////////
	      
	      if(aLayer->IsX()) 
		{
		  // Fill the clusters
		  for (int j=0; j<ClusterLayerSize;j++)
		    {
		      XClusterXlayer.push_back(ClusterLayer[j]);
		      ZClusterXlayer.push_back(height);
		    }
		} 
	      else
		{ 
		  // Fill the clusters
		  for (int j=0; j<ClusterLayerSize;j++)
		    {
		      XClusterYlayer.push_back(ClusterLayer[j]);
		      ZClusterYlayer.push_back(height);
		    }
		}
	    }
	}
      int NumClusX = XClusterXlayer.size();
      int NumClusY = XClusterYlayer.size();  
      
      XClusters->Set(NumClusX);
      YClusters->Set(NumClusY);
      
      for (int i =0;i<NumClusX;i++)
	{
	  XClusters->SetPoint(i,XClusterXlayer[i],ZClusterXlayer[i]);
	}
      for (int i =0;i<NumClusY;i++)
	{
	  YClusters->SetPoint(i,XClusterYlayer[i],ZClusterYlayer[i]);
	}
      //////////////////////////////////////////////////
      gDirectory->Delete("fun");
      TF1 *fun = new TF1("fun","pol1",0,40);

      double Xextrapolated, Yextrapolated;
      double Ax = 0;
      double Bx = 0;      
      double Ay = 0;
      double By = 0;
      
      if(NumClusX > 1 && NumClusY > 1)
	{
	  XClusters->Fit("fun","QR");
	  Ax = fun->GetParameter(0);
	  Bx = fun->GetParameter(1);
	  YClusters->Fit("fun","QR");
	  Ay = fun->GetParameter(0);
	  By = fun->GetParameter(1);
	  
	  ti.Reset();
	  
	  while ( ( key = (TObjString*)ti.Next() ) ) 
	    {
	      
	      Layer *aLayer = ((Layer*) myGeometry->GetValue(key));
	      if(LV!="All" && key->GetString()!=LV) continue;
	      //	      if(LV=="Y6") continue;
	      double Zlayer = aLayer->GetHeight();
	      int LayerNumHits = myEvent->GetLayerNumHits(key->GetString());
	      
	      Yextrapolated = (Zlayer - Ay)/By;
	      Xextrapolated = (Zlayer - Ax)/Bx;
	      bool ActiveRegion=false;
	      if((aLayer->IsX() && (aLayer->checkActiveArea(Xextrapolated,Yextrapolated,1.0))) ||
		 (aLayer->IsY() && (aLayer->checkActiveArea(Yextrapolated,Xextrapolated,1.0))))		    
		{ 
		  ActiveRegion=true;
		  XInActiveArea.push_back(Xextrapolated);
		  YInActiveArea.push_back(Yextrapolated);
		  aLayer->AddHitInActiveArea();
		}
 	      else
		{
		  XNotInActiveArea.push_back(Xextrapolated);
		  YNotInActiveArea.push_back(Yextrapolated);
		}
	      
	      if(ActiveRegion && LayerNumHits==0)
		{
		  XMissedHits.push_back(Xextrapolated);
		  YMissedHits.push_back(Yextrapolated);
		  aLayer->AddMissedHit();
		}
	      else if(ActiveRegion && LayerNumHits>0)
		{
		  std::vector<double> ClusterLayer = myEvent->GetClusters(aLayer->GetLayerName());
		  int NumberOfClusterPerLayer      = ClusterLayer.size();
		  
		  if(aLayer->IsX() && NumberOfClusterPerLayer==1) 
		    XResiduals->Fill(Xextrapolated-ClusterLayer[0]);
		  else if(aLayer->IsY() && NumberOfClusterPerLayer==1)
		    YResiduals->Fill(Yextrapolated-ClusterLayer[0]);
		}
	    }
	}
    }
  
  int NNInAA = XNotInActiveArea.size();
  int NInAA  = XInActiveArea.size();
  int NMH    = XMissedHits.size();

  double *x1 = new double[NNInAA];
  double *x2 = new double[NInAA];
  double *x3 = new double[NMH];
  double *y1 = new double[NNInAA];
  double *y2 = new double[NInAA];
  double *y3 = new double[NMH];

  for(int i =0; i<NNInAA; i++)
    {
      x1[i] = XNotInActiveArea[i];
      y1[i] = YNotInActiveArea[i];
    }
  
  for(int i =0; i<NInAA; i++)
    {
      x2[i] = XInActiveArea[i];
      y2[i] = YInActiveArea[i];
    }
 
  for(int i =0; i < NMH; i++)
    {
      x3[i] = XMissedHits[i];
      y3[i] = YMissedHits[i];
    }

  double FakeCoordinate[2]={0,40};
  TGraph *FakeAxis = new TGraph(2,FakeCoordinate,FakeCoordinate);
  TGraph *XYRegionNAA = new TGraph(NNInAA,x1,y1);
  TGraph *XYRegionAA  = new TGraph(NInAA,x2,y2);
  TGraph *XYRegionMH  = new TGraph(NMH,x3,y3);
  
	
  TCanvas *c1;
  c1 = new TCanvas("c1","c1",500,500);
  XYRegionAA->SetMarkerColor(2);
  XYRegionNAA->SetMarkerColor(3);
  XYRegionMH->SetMarkerStyle(28);

  FakeAxis->Draw("ap");
  XYRegionAA->Draw("p");
  XYRegionNAA->Draw("p");
  XYRegionMH->Draw("p");

  ti.Reset();
  Layer *aLayer;
  std::cout<<" LAYER NAME | EFFICIENCY | INEFFICIENCY"<<std::endl;
  while ( ( key = (TObjString*)ti.Next() ) ) 
    {
      aLayer = ((Layer*) myGeometry->GetValue(key));
      std::cout<<" |   "<<aLayer->GetLayerName()<<"    |     ";
      std::cout << std::setw(5) << std::setprecision(2) << aLayer->GetEfficiency()*100.0 << " %  |"
                << std::setw(6) << "   % " << aLayer->GetInefficiency()*100.0;
      //      printf("  %3.2f \%  |   %3.2f   \% ", (float) aLayer->GetEfficiency()*100.0,(float) aLayer->GetInefficiency()*100.0);
      std::cout<<"   | "<<std::endl;
    }

  TCanvas *c2 = new TCanvas("c2","c2",500,500);
  c2->Divide(2,1);
  c2->cd(1);
  XResiduals->Draw();
  c2->cd(2);
  YResiduals->Draw();
}

