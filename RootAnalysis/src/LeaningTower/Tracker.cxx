#include "Tracker.h"
#include "Layer.h"
#include "TObjString.h"
#include "TCanvas.h"
#include "TGraph.h"

#include <iomanip>
#include <iostream>
#include <fstream>
#include <string>

Tracker::Tracker()
{
  myGeometry = new TMap();
  TOWER=false;
}


void Tracker::loadGeometry(TString filename) {
  if ( !myGeometry->IsEmpty() ) {
    std::cout << "WARNING: Geometry not empty!  Maybe you filled it before?" << std::endl;
    return;
  }
  std::ifstream fin(filename);
  if ( !fin ) {
      std::cout << "File " << filename << " could not be opened!" << std::endl;
      return;
  }
  TString layer;
  double value;
  while ( 1 ) 
    {
      fin >> layer >> value;
      if ( !fin.good() ) break;
      Layer *aLayer = new Layer();
      aLayer->SetLayerName(layer);
      aLayer->SetHeight(value);
      myGeometry->Add(new TObjString(layer),aLayer);
    }
  fin.close();
  myGeometry->Print();
}

void Tracker::Display()
{
  TCanvas *EventDisplayC = (TCanvas*) gROOT->FindObject("EventDisplayC");
  
  TMapIter ti(myGeometry);
  TObjString* key;
  
  double x[2] = { -4, 40};
  double y[2]; 
  if(TOWER)
    {
      y[0]=0;
      y[1]=70;
    }
  else
    {
      y[0]=-110;
      y[1]=110;
    }
  //  double y[2] = { -0, 650}; //TOWER
  //  double y[2] = { -110, 100}; //STACK

  TGraph *Border = new TGraph(2, x, y);
  Border->GetXaxis()->SetTitle("position [cm]");
  Border->GetYaxis()->SetTitle("position in stack/cm");

  EventDisplayC->Divide(2,1);
  EventDisplayC->cd(1);
  Border->Draw("AP");
  EventDisplayC->cd(2);
  Border->Draw("AP");
  
  while ( ( key = (TObjString*)ti.Next() ) ) 
    {
      Layer *aLayer = ((Layer*) myGeometry->GetValue(key));
      
      
      if(aLayer->IsX()) 
	{
	  EventDisplayC->cd(1);
	  aLayer->DrawLayer();	      
	  EventDisplayC->cd(2);
	  aLayer->DrawGhostLayer();	
	}
      else
	{
	  EventDisplayC->cd(2);
	  aLayer->DrawLayer();	      
	  EventDisplayC->cd(1);
	  aLayer->DrawGhostLayer();
	}

      
      //      std::cout<<((Layer*) myGeometry->GetValue(key))->GetLayerName()<<std::endl;
    }
  if(!TOWER)
    {
      TLine *ScintillatorT  = new TLine(-2, 89.5,38,89.5);
      TLine *ScintillatorM  = new TLine(-2, 0.0,38,0.0);
      TLine *ScintillatorB  = new TLine(-2, -103.1,38,-103.1);
      
      TLine *ScintillatorTL = new TLine(-2, 89.5, 17.75,89.5);
      TLine *ScintillatorTR = new TLine(18.25, 89.5, 38,89.5); 
      TLine *ScintillatorML = new TLine(-2, 0.0, 17.75,0.0);
      TLine *ScintillatorMR = new TLine(18.25, 0.0, 38,0.0);
      TLine *ScintillatorBL = new TLine(-2, -103.1, 17.75,-103.1);
      TLine *ScintillatorBR = new TLine(18.25, -103.1, 38,-103.1);
      
      ScintillatorTL->SetLineColor(4);
      ScintillatorTR->SetLineColor(4);
      ScintillatorML->SetLineColor(4);
      ScintillatorMR->SetLineColor(4);
      ScintillatorBL->SetLineColor(4);
      ScintillatorBR->SetLineColor(4);  
      ScintillatorT->SetLineColor(4);  
      ScintillatorM->SetLineColor(4);
      ScintillatorB->SetLineColor(4); 
      
      EventDisplayC->cd(1);
      ScintillatorTL->Draw();
      ScintillatorTR->Draw();
      ScintillatorML->Draw();
      ScintillatorMR->Draw();
      ScintillatorBL->Draw();
      ScintillatorBR->Draw();
      
      EventDisplayC->cd(2);
      ScintillatorT->Draw();
      ScintillatorM->Draw();
      ScintillatorB->Draw();
    }
}
