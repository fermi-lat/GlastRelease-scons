#include <iomanip>
#include <iostream>
#include <fstream>
#include <string>

#include "TGraph.h"
#include "TTree.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TString.h"

TMap* myGeometry = new TMap();

void EventDisplay(TString="MyRootFile.root", int=0);
void loadGeometry(TString filename="src/LeaningTower/geometry/stack2geometry.txt");
void getGeometry(TString layer);

double getGeometry(TString layer) {
    if ( myGeometry->IsEmpty() ) {
        std::cerr << "ERROR:  Geometry is empty!  Maybe you forgot to load it?" << std::endl;
        return;
    }
    if ( !myGeometry->FindObject(layer) ) {
        std::cerr << "ERROR: Layer " << layer << " not found in geometry map!" << std::endl;
        return 1E15;
    }
    return atof(((TObjString*)myGeometry->GetValue(layer))->GetName());
}

void loadGeometry(TString filename) {
    if ( !myGeometry->IsEmpty() ) {
        std::cout << "WARNING: Geometry not empty!  Maybe you filled it before?" << std::endl;
        return;
    }
    std::ifstream fin(filename, ios::in);
    TString layer;
    TString value;
    while ( 1 ) {
        fin >> layer >> value;
        if ( !fin.good() )
            break;
        myGeometry->Add(new TObjString(layer), new TObjString(value));
    }
    fin.close();
    myGeometry->Print();
}

void EventDisplay(TString filename, int firstEvent) {
    myFile = new TFile(filename);
    myTree = (TTree*)myFile->Get("Header");
    int numEvents = myTree->GetEntries();
    std::cout << "Number of Events: " << numEvents << std::endl;
    TCanvas evDisp("c1","Event Display", 700, 800);
    Int_t EventId;
    Int_t RunId;
    Int_t TkrTotalNumHits;  
    myTree->SetBranchAddress("EventId",&EventId);
    myTree->SetBranchAddress("RunId",&RunId);
    myTree->SetBranchAddress("TkrTotalNumHits",&TkrTotalNumHits);

    for ( int iEvent=firstEvent; iEvent<numEvents; ++iEvent ) {
        myTree->GetEntry(iEvent);
    
        std::cout << "EventId: " << EventId
                  << " RunId: " << RunId
                  << " TkrTotalNumHits: "  << TkrTotalNumHits
                  << std::endl;
        if ( myGeometry->IsEmpty() ) {
            std::cerr << "ERROR:  Geometry is empty!"
                      << "  Maybe you forgot to load it?"
                      << std::endl;
            return;
        }
        TMapIter ti(myGeometry);
        TObjString* key;
        TTree* layer;
        double x[4] = { 0, 1550, 1550, 0 };
        double y[4] = { -65, -65, 65, 65 };
        TGraph* tg = new TGraph(4, x, y);
        tg->GetXaxis()->SetTitle("strip Id");
        tg->GetYaxis()->SetTitle("position in stack/cm");
        tg->Draw("AP");
        TLine* l;
        TText* t = new TText(0, 70, TString("EventId: ") + EventId + " TkrTotalNumHits: " + TkrTotalNumHits);
        t->SetTextAlign();
        t->SetTextSize(0.04);
        t->Draw();
        while ( key = (TObjString*)ti.Next() ) {
            TString VL = key->String();
            double pos = getGeometry(VL);
            l = new TLine(0, pos, 1535, pos);
            l->Draw();
            layer = (TTree*)myFile->Get(TString("Layer")+VL);
            Int_t TkrNumHits;
            Int_t TkrHits[128];
            layer->SetBranchAddress("TkrNumHits",&TkrNumHits);
            layer->SetBranchAddress("TkrHits",TkrHits);
            layer->GetEntry(iEvent);
            t = new TText(1550, pos, VL + " (" + TkrNumHits + ")" );
            t->SetTextAlign(02);
            t->SetTextSize(0.02);
            t->Draw();
            if ( TkrNumHits > 0 ) {
                Double_t xpos[128];
                Double_t zpos[128];
                for ( int i=0; i<TkrNumHits; ++i ) {
                    xpos[i] = TkrHits[i];
                    zpos[i] = pos;
                }
                tg = new TGraph(TkrNumHits, xpos, zpos);
                tg->Draw("*");
            }
        }
        evDisp.Update();
        std::cout << "Hit return to continue (q to quit): ";
        char c = getchar();
        if ( c == 'q' )
            return;
    }
}
