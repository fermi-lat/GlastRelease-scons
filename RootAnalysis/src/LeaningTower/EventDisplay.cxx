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
TText TriggerReqText[36][2];

TGraph *XTrack;
TGraph *YTrack;

TGraph *XClusters;
TGraph *YClusters;
TLine anXtrack,anYtrack;

int EventId, RunId,TkrTotalNumHits;
TH1D *CHI2;
static const double MaxChi2 = 1;

//////////////////////////////////////////////////
#define TOWER 1


TLine *Track(TGraph *XY)
{
  int N = XY->GetN();
  double *X = XY->GetX();
  double *Y = XY->GetY();
  double Xm = 0;
  double Ym = 0;
  for (int i = 0 ; i<N; i++)
    {
      Xm+=X[i];
      Ym+=Y[i];
    }
  Xm/=N;
  Ym/=N;
  double Xr,Yr;
  double A=0.0;
  double B=0.0;
  for (int i = 0 ; i<N; i++)
    {
      Xr = (Xm-X[i]);
      Yr = (Ym-Y[i]);
      A+=(Xr*Yr);
      B+=pow(Yr,2.0)-pow(Xr,2.0);
    }
  double Theta = -0.5*TMath::ATan2(2.0*A,B);

  TLine *track = new TLine(Xm,Ym,Xm+10* TMath::Cos(Theta),Ym+10* TMath::Sin(Theta));
  return track;
  
}

TLine Recon(TGraph *XY)
{

  double Y1=0.0;
  double Y2=70.0;
  if(TOWER)
    {
      Y1=0.0;
      Y2=70.0;
    }
  int N = XY->GetN();
  if ( N < 2 )    
      return TLine(0,0,0,0);

  double *X = XY->GetX();
  double *Y = XY->GetY();
  double *newX = new double[N-1];
  double *newY = new double[N-1];
  double MinChi2 = 1e6;
  double A=0;
  double B=0;
  TLine bestFit;

  if(N<=3)
    {
      TGraph g(N,Y,X);
      g.Fit("pol1","Q");
      MinChi2 = ((TF1*) g.FindObject("pol1"))->GetChisquare();      
      A = ((TF1*) g.FindObject("pol1"))->GetParameter(0);
      B = ((TF1*) g.FindObject("pol1"))->GetParameter(1);
    }
  else
    {
      for (int j = 0; j<N ; j++)
	{
	  int l=0;
	  for(int i=0;i<N;i++)
	    {
	      if(i!=j) 
		{
		  
		  newX[l]  = X[i];
		  newY[l++]= Y[i];
		}
	    }
	  TGraph g(N-1,newY,newX);
	  g.Fit("pol1","Q");
	  double chi2 = ((TF1*) g.FindObject("pol1"))->GetChisquare();
	  //	  std::cout<<chi2<<" "<<MinChi2<<std::endl;
	  if(chi2<MinChi2)
	    {
	      MinChi2=chi2;
	      A = ((TF1*) g.FindObject("pol1"))->GetParameter(0);
	      B = ((TF1*) g.FindObject("pol1"))->GetParameter(1);
	    }
	}
    }
  std::cout<<"chisqr " <<MinChi2<<std::endl;
  if(MinChi2>MaxChi2)
    return TLine(0,0,0,0);
  

  double x1,x2,y1,y2;
  x1 = A+B*Y1;//0.0;
  y1 = Y1;//A;
  x2 = A+B*Y2;//W;
  y2 = Y2;//;
  //  std::cout<<x1<<" "<<y1<<" "<<x2<<" "<<y2<<std::endl;
  TLine track(x1,y1,x2,y2);
  track.SetLineColor(2);
  return track;
  
}

//////////////////////////////////////////////////
bool MeritFactor(TLine atrack)
{
  double x1 = atrack.GetX1();
  double x2 = atrack.GetX2();
  double y1 = atrack.GetY1();
  double y2 = atrack.GetY2();
  bool myFlag=true;
  if((x1 == 0.0 && x2 == 0.0) &&(y1 == 0.0 && y2 ==0.0)) 
    myFlag=false;
  
  return  myFlag;
}


//////////////////////////////////////////////////
double ExtrapolateCoordinate(TLine track, double Z)
{
  double x1 = track.GetX1();
  double x2 = track.GetX2();
  double y1 = track.GetY1();
  double y2 = track.GetY2();
  
  double XX = (x2==x1)? x1 : x1 + (Z-y1)/(y2-y1)*(x2-x1);
  return XX;
}


//////////////////////////////////////////////////
void InitializeED(TString filename = "MyRootFile.root")
{
  //  Initialize(filename);
  gStyle->SetCanvasColor(10);
  tracker = new Tracker();

  if(TOWER)
    {  
      tracker->loadGeometry(gSystem->ExpandPathName("$ROOTANALYSISROOT/src/LeaningTower/geometry/Tower0Geometry.txt"));
      tracker->IsTower(TOWER);
    }
  else
    tracker->loadGeometry(gSystem->ExpandPathName("$ROOTANALYSISROOT/src/LeaningTower/geometry/stack2geometry.txt"));
  
  myEvent = new Event(filename,(TMap *) tracker->GetGeometry());
  
  EventDisplayC = new TCanvas("EventDisplayC"  ,"EventDisplayC",1100,800);  
  tracker->Display();
  
  for ( int i=0; i<36; ++i ) {
    TriggerReqText[i][0].SetTextAlign(32);
    TriggerReqText[i][1].SetTextAlign(02);
    LabelNumHits[i].SetTextAlign(02);
    TriggerReqText[i][0].SetTextSize(0.02);
    TriggerReqText[i][1].SetTextSize(0.02);
    LabelNumHits[i].SetTextSize(0.02);
  }

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

void DisplayEvent(int NumEvent=-1)
{
  static int realEvent = -1;
  if ( NumEvent < 0 )
      ++realEvent;
  else
      realEvent = NumEvent;

  myEvent->Go(realEvent);
  
  int TkrTotalNumHits = myEvent->GetTkrTotalNumHits();


  std::cout<<" TkrTotalNumHits = "<<TkrTotalNumHits<<std::endl;
  std::cout<<" EventId         = "<<myEvent->GetEventId()<<std::endl;
  std::cout<<" RunId           = "<<myEvent->GetRunId()<<std::endl;
  std::cout<<" EbfTime         = "<<myEvent->GetEbfTime()<<std::endl;

  
  TMap *myGeometry = tracker->GetGeometry();
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
  
  TText base;
  base.SetTextAlign(02);
  base.SetTextSize(0.02);

  // let's do Recon here



  TMapIter ti(myGeometry);
  TObjString* key;
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

      if(aLayer->IsX()) 
	EventDisplayC->cd(1);      
      else
	EventDisplayC->cd(2);

      TriggerReqText[i][0].SetText(-0.8, height, aLayer->GetTriggerReq(false) ? "x" : ".");
      TriggerReqText[i][1].SetText(36, height, aLayer->GetTriggerReq(true) ? "x" : ".");
      LabelNumHits[i].SetText(41.0, height, title);
      TriggerReqText[i][0].Draw();
      TriggerReqText[i][1].Draw();
      LabelNumHits[i].Draw();
      i++;      

      if(LayerNumHits>0)
	{
	  //////////////////////////////////////////////////
	  //                   Clustering
	  std::vector<double> ClusterLayer = myEvent->GetClusters(LayerName);
	  int ClusterLayerSize = ClusterLayer.size();
	  //////////////////////////////////////////////////
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
  // Display the hits
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
  // RECON
  
  EventDisplayC->cd(1);
  if ( NumClusX > 0 )
      XClusters->Draw("P");
  anXtrack = Recon(XClusters);
  anXtrack.Draw();
  
  EventDisplayC->cd(2);
  if ( NumClusY > 0 )
      YClusters->Draw("P");
  anYtrack = Recon(YClusters);
  anYtrack.Draw();
  
  EventDisplayC->cd();
  EventDisplayC->Update();
  
}


void IneffAnalysis(int LastEvent, TString LV="All")
{

  //  CHI2 = new TH1D("CHI2","CHI2",100,0,MaxChi2);

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
  int SelectedEvents=0;
  int TotalEvents=0;
  for(int ev = 1;ev <= LastEvent; ev++)
    { 
      TotalEvents++;
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
    
      double Xextrapolated, Yextrapolated;
      
      if(NumClusX > 1 && NumClusY > 1)
	{
	  TLine Xtrack = Recon(XClusters);
	  TLine Ytrack = Recon(YClusters);
	  
	  if(MeritFactor(Xtrack) && MeritFactor(Ytrack))
	    {
	      ti.Reset();
	      SelectedEvents++;
	      
	      while ( ( key = (TObjString*)ti.Next() ) ) 
		{
		  Layer *aLayer = ((Layer*) myGeometry->GetValue(key));
		  if(LV!="All" && key->GetString()!=LV) continue;
		  //	      if(LV=="Y6") continue;
		  double Zlayer = aLayer->GetHeight();
		  int LayerNumHits = myEvent->GetLayerNumHits(key->GetString());
		  
		  Xextrapolated = ExtrapolateCoordinate(Xtrack,Zlayer);
		  Yextrapolated = ExtrapolateCoordinate(Ytrack,Zlayer);
		  
		  //		  std::cout<<"Z = "<<Zlayer<<" Xextrapolated = "<<Xextrapolated<<" Yextrapolated = "<<Yextrapolated<<std::endl;
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
      std::cout << std::setw(6) << std::setprecision(3) << aLayer->GetEfficiency()*100.0 << " %  |"
                << std::setw(6) << "   % " << aLayer->GetInefficiency()*100.0;
      //      printf("  %3.2f \%  |   %3.2f   \% ", (float) aLayer->GetEfficiency()*100.0,(float) aLayer->GetInefficiency()*100.0);
      std::cout<<"   | "<<std::endl;
    }
  std::cout<<" SELECTED EVENTS = "<<SelectedEvents<<" / "<<TotalEvents<<std::endl;
  TCanvas *c2 = new TCanvas("c2","c2",500,500);
  c2->Divide(2,1);
  c2->cd(1);
  XResiduals->Draw();
  c2->cd(2);
  YResiduals->Draw();
  c2->cd();
  /*
    TCanvas *c3 = new TCanvas();
    CHI2->Draw();
  */
}

