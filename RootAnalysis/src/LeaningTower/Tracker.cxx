#include "Tracker.h"

#include "Layer.h"

#include "TObjString.h"
#include "TGraph.h"

#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

Tracker::Tracker()
{
  myGeometry = new TMap();
  TOWER=false;
}

Tracker::~Tracker() {
    myGeometry->DeleteAll();  // deletes the keys and values of the map
    delete myGeometry;        // deletes the map
}


void Tracker::loadGeometry(TString filename) {
  if ( !myGeometry->IsEmpty() ) {
    std::cout << "WARNING: Geometry not empty!  Maybe you filled it before?"
              << std::endl;
    return;
  }
  std::ifstream fin(filename);
  if ( !fin ) {
      std::cout << "File " << filename << " could not be opened!" << std::endl;
      return;
  }
  std::string line;
  while ( !fin.eof() ) {
      std::getline(fin, line);
      TString layer;
      double z, y, x;
      z = y = x = 0;
      std::istringstream ist(line);
      ist >> layer >> z >> y >> x;
      if ( layer == "" )
          continue;
      Layer* aLayer = new Layer(layer, z, y, x);
      myGeometry->Add(new TObjString(layer), aLayer);
      std::cout << "layer " << layer << std::endl;
  }

  fin.close();
  myGeometry->Print();
}

void Tracker::loadFitting(TString filename) {
  if ( myGeometry->IsEmpty() ) {
    std::cout << "WARNING: Geometry is empty!" << std::endl;
    return;
  }
  std::ifstream fin(filename);
  if ( !fin ) {
      std::cout << "File " << filename << " could not be opened!" << std::endl;
      return;
  }
  TString name;
  double value;
  while ( 1 ) {
      fin >> name >> value;
      if ( !fin.good() ) break;
      std::vector<TString> v;
      for ( int i=0; i<value; ++i ) {
          TString dummy;
          fin >> dummy;
          v.push_back(dummy);
      }
      Layer* l = (Layer*)myGeometry->GetValue(name);
      l->SetPlanesForFittingCol(v);
    }
  fin.close();
}

std::vector<TString> Tracker::GetPlaneNameCol(const TString view,
                                              const bool onlyTheGood) const {
    if ( view == "X" )
        return GetPlaneNameCol(0, onlyTheGood);
    else if ( view == "Y" )
        return GetPlaneNameCol(1, onlyTheGood);
    else
        std::cerr << "Tracker::GetPlaneNameCol: view = " << view << std::endl;

    std::exit(42);
    return std::vector<TString>();
}

std::vector<TString> Tracker::GetPlaneNameCol(const int view,
                                              const bool onlyTheGood) const {
    TMapIter ti(myGeometry);
    TObjString* key;
    std::vector<TString> v;
    while ( ( key = (TObjString*)ti.Next() ) ) {
        Layer* l = (Layer*)myGeometry->GetValue(key);
        if ( view != l->GetView() )
            continue;
        if ( onlyTheGood && l->GetPlanesForFittingCol().size() == 0 )
            continue;
        v.push_back(l->GetLayerName());
    }
    return v;
}

void Tracker::Display(TCanvas* ed) {
  TMapIter ti(myGeometry);
  TObjString* key;
  
  double x[2] = { -40, 400};
  double y[2]; 
  if(TOWER)
    {
      y[0]=0;
      y[1]=700;
    }
  else
    {
      y[0]=-1100;
      y[1]=1100;
    }
  //  double y[2] = { -0, 6500}; //TOWER
  //  double y[2] = { -1100, 1000}; //STACK

  TGraph *Border = new TGraph(2, x, y);
  Border->GetXaxis()->SetTitle("position/mm");
  Border->GetYaxis()->SetTitle("position in stack/mm");

  ed->Divide(2,1);
  ed->cd(1);
  Border->Draw("AP");
  ed->cd(2);
  Border->Draw("AP");
  
  while ( ( key = (TObjString*)ti.Next() ) ) 
    {
      Layer *aLayer = ((Layer*) myGeometry->GetValue(key));
      
      
      if(aLayer->IsX()) 
	{
	  ed->cd(1);
	  aLayer->DrawLayer();	      
	  ed->cd(2);
	  aLayer->DrawGhostLayer();	
	}
      else
	{
	  ed->cd(2);
	  aLayer->DrawLayer();	      
	  ed->cd(1);
	  aLayer->DrawGhostLayer();
	}
    }
  if(!TOWER)
    {
      TLine *ScintillatorT  = new TLine(-20,    895,  380,     895);
      TLine *ScintillatorM  = new TLine(-20,      0,  380,       0);
      TLine *ScintillatorB  = new TLine(-20,  -1031,  380,   -1031);
      
      TLine *ScintillatorTL = new TLine(-20,    895,  177.5,   895);
      TLine *ScintillatorTR = new TLine(182.5,  895,  380,     895); 
      TLine *ScintillatorML = new TLine(-20,      0,  177.5,     0);
      TLine *ScintillatorMR = new TLine(182.5,    0,  380,       0);
      TLine *ScintillatorBL = new TLine(-20,   -1031, 177.5, -1031);
      TLine *ScintillatorBR = new TLine(182.5, -1031, 380,   -1031);
      
      ScintillatorTL->SetLineColor(4);
      ScintillatorTR->SetLineColor(4);
      ScintillatorML->SetLineColor(4);
      ScintillatorMR->SetLineColor(4);
      ScintillatorBL->SetLineColor(4);
      ScintillatorBR->SetLineColor(4);  
      ScintillatorT->SetLineColor(4);  
      ScintillatorM->SetLineColor(4);
      ScintillatorB->SetLineColor(4); 
      
      ed->cd(1);
      ScintillatorTL->Draw();
      ScintillatorTR->Draw();
      ScintillatorML->Draw();
      ScintillatorMR->Draw();
      ScintillatorBL->Draw();
      ScintillatorBR->Draw();
      
      ed->cd(2);
      ScintillatorT->Draw();
      ScintillatorM->Draw();
      ScintillatorB->Draw();
    }
  ed->cd();
}

ClassImp(Tracker)
