#include "Tracker.h"

#include "Layer.h"

#include "TObjString.h"
#include "TGraph.h"

#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

Tracker::Tracker() {
    // myGeometry was a TMap.  The key was name, and the object was a Layer
    // object.  A TMap is stored in arbitrary order, and cannot be sorted.
    // There is also no other advantage in having a TMap, as the information is
    // redundant (the name is contained in the Layer).

    // Thus, let's use a TList.  A TList can be filled ordered, and moreover can
    // be sorted, if neccessary.

    // Actually, it could also be a std::vector or std::list.  However, to
    // retrieve a particular entry one would have to define a predicate.  Would
    // root compile that?

    myGeometry = new TList();
    TOWER=false;
}

Tracker::~Tracker() {
    // Remove all objects from the list AND delete all heap based objects.
    myGeometry->Delete();
    delete myGeometry;        // deletes the TList
}


void Tracker::loadGeometry(TString filename) {
    gSystem->ExpandPathName(filename);
    if ( !myGeometry->IsEmpty() ) {
        std::cout << "WARNING: Geometry not empty!  Maybe you filled it before?"
                  << std::endl;
        return;
    }
    std::ifstream fin(filename);
    if ( !fin ) {
        std::cout << "File " << filename << " could not be opened!" <<std::endl;
        return;
    }
    std::string line;
    while ( !fin.eof() ) {
        std::getline(fin, line);
        TString layer;
        double z, y, x, rz, ry, rx;
        z = y = x = rz = ry = rx = 0;
        std::istringstream ist(line);
        ist >> layer >> z >> y >> x >> rz >> ry >> rx;
        if ( layer == "" )
            continue;
        Layer* aLayer = new Layer(layer, z, y, x, rz, ry, rx);
        myGeometry->AddLast(aLayer);
    }

    fin.close();
    myGeometry->Print();
}

void Tracker::loadFitting(TString filename) {
    gSystem->ExpandPathName(filename);
    if ( myGeometry->IsEmpty() ) {
        std::cout << "WARNING: Geometry is empty!" << std::endl;
        return;
    }
    std::ifstream fin(filename);
    if ( !fin ) {
        std::cout << "File " << filename << " could not be opened!" <<std::endl;
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
        Layer* l = (Layer*)myGeometry->FindObject(name);
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
    TIter next(myGeometry);
    std::vector<TString> v;
    while ( Layer* l = (Layer*)next() ) {
        if ( view != l->GetView() )
            continue;
        if ( onlyTheGood && l->GetPlanesForFittingCol().size() == 0 )
            continue;
        v.push_back(l->GetName());
    }
    return v;
}

void Tracker::Display(TCanvas* ed) {
    double x[2] = { -40, 400};
    double y[2]; 
    if ( TOWER ) {
        y[0] = 0;
        y[1] = 700;
    }
    else {
        y[0] = -1100;
        y[1] = 1100;
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
  
    TIter next(myGeometry);
    while ( Layer* aLayer = (Layer*)next() ) {
        if ( aLayer->IsX() ) {
            ed->cd(1);
            aLayer->DrawLayer();	      
            ed->cd(2);
            aLayer->DrawGhostLayer();	
	}
        else {
            ed->cd(2);
            aLayer->DrawLayer();	      
            ed->cd(1);
            aLayer->DrawGhostLayer();
	}
    }
    if ( !TOWER ) {
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
