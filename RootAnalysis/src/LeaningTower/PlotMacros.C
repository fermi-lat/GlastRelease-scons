#include <iomanip>
#include <iostream>

gROOT->Reset();
TFile myFile;
TFile myOutFile;
TTree *myTree;
TString TID;

const int NUM_ROC = 2;
const int NUM_LAYER = 18;
const int NUM_VIEW = 2;
void help();
void Init(char *filename="MyRootFile.root", char* f="MyAnalysisFile.root")
{
  gStyle->SetCanvasColor(10);
  gStyle->SetFillColor(10);
  
  TFile myFile(filename);
  myTree = (TTree*) myFile.Get("myTree");
  myTree->StartViewer();
  std::cout<<TrayId(1,0)<<std::endl;
  TFile myOutFile(f, "recreate");
  help();
}

void help()
{
  std::cout<<" ******************** Analysis Program for Tray data ********************"<<std::endl;
  std::cout<<" ---------- DIGI ANALYSIS:     "<<std::endl;
  std::cout<<" TrayId(l , v) "<<std::endl;
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
  
  TString ts="[";
  //  TString *trayid;
  ts +=l ;
  ts+="][";
  ts+=v;
  ts+="]";
  TID=ts;
  return ts;
}


void PlotToT_0(TString myCuts="")
{
  TString var="ToT[0]";
  var+=TID;

  gDirectory->Delete("ToT0");
  TH1D *ToT0 = new TH1D("ToT0",var,255,0,255);
  var+=">>ToT0";
  myTree->Draw(var,myCuts);
}

void PlotToT_1(TString myCuts="")
{
  TString var="ToT[1]";
  var+=TID;

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
  var+=TID;
  gDirectory->Delete("HitMap");
  TH1D *HitMap = new TH1D("HitMap",var,1536,0,1536);
  HitMap->SetLineColor(4);
  var+=">>HitMap";
  myTree->Draw(var,myCuts);
  HitMap->SetEntries(HitMap->GetEntries()/128);
  std::cout<<" Number of selected Events: "<<HitMap->GetEntries()<<std::endl;
}

void PlotNumHitsLayer(TString myCuts="")
{
  gDirectory->Delete("NumHitsCan");
  TCanvas *NumHitsCan = new TCanvas("NumHitsCan","Num Hits",500,500);
  TString var="TkrNumHits";
  var+=TID;
  gDirectory->Delete("TkrNumHits");
  //  TH1D *TkrNumHits = new TH1D("TkrNumHits",var,1536,0,1536);
  var+=">>TkrNumHits";
  myTree->Draw(var,myCuts);
  TkrNumHits->SetFillColor(3);

}

void FindEvents(TString myCuts="")
{
  myTree->Scan("EventId",myCuts);
}

void PlotTotTot(TString myCuts="")
{ 
  gDirectory->Delete("ToTToTCan");
  TCanvas *ToTToTCan = new TCanvas("ToTToTCan","ToT vs ToT",500,500);
  TString var0="ToT[0]";
  TString var1="ToT[1]";
  var0+=TID;
  var1+=TID;
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

bool DEBUG = 0;
void debug() {
    DEBUG = !DEBUG;
}

void Analyze() {
    // Purpose and Method:  Analyze the TTree for oddities

    myOutFile.cd();

    TH1F* TkrDigiHitMap[NUM_LAYER][NUM_VIEW];
    for ( int l=0; l<NUM_LAYER; ++l )
        for ( int v=0; v<NUM_VIEW; ++v ) {
            TString s = TrayId(l, v);
            TkrDigiHitMap[l][v] = new TH1F("TkrDigiHitMap"+s,
                                           "TkrDigi hit map"+s+";strip id",
                                           1536, -0.5, 1535.5);
        }

    int nEntries = (int)myTree->GetEntries();
    std::cout << "myTree entries: " << nEntries << std::endl;

    Int_t EventId;
    Int_t RunId;
    Int_t ToT[NUM_ROC][NUM_LAYER][NUM_VIEW];
    Int_t TkrTotalNumHits;
    Int_t TkrNumHits[NUM_LAYER][NUM_VIEW];
    Int_t TkrHits[NUM_LAYER][NUM_VIEW][128];

    myTree->SetBranchAddress("EventId", &EventId);
    myTree->SetBranchAddress("RunId", &RunId);
    myTree->SetBranchAddress("ToT", ToT);
    myTree->SetBranchAddress("TkrTotalNumHits", &TkrTotalNumHits);
    myTree->SetBranchAddress("TkrNumHits", TkrNumHits);
    myTree->SetBranchAddress("TkrHits", TkrHits);

    Int_t TkrNumClus;
    Int_t TkrClusPlane[128];
    Int_t TkrClusView[128];
    Double_t TkrClusX[128];
    Double_t TkrClusY[128];
    Double_t TkrClusZ[128];
    Int_t TkrNumTracks;
    Int_t TkrTrk1Clusters[128];
    myTree->SetBranchAddress("TkrNumClus",&TkrNumClus);
    myTree->SetBranchAddress("TkrClusX",TkrClusX);
    myTree->SetBranchAddress("TkrClusY",TkrClusY);
    myTree->SetBranchAddress("TkrClusZ",TkrClusZ);
    myTree->SetBranchAddress("TkrClusPlane",TkrClusPlane);
    myTree->SetBranchAddress("TkrClusView",TkrClusView);
    myTree->SetBranchAddress("TkrNumTracks",&TkrNumTracks);
    myTree->SetBranchAddress("TkrTrk1Clusters",TkrTrk1Clusters);

    //    nEntries = 1;
    for ( int iEvent=0; iEvent<nEntries; ++iEvent ) {
        myTree->GetEntry(iEvent);

        if ( DEBUG ) {
            std::cout << "EventId " << EventId;
            std::cout << " RunId " << RunId << std::endl;
            std::cout << "TkrTotalNumHits " << TkrTotalNumHits << std::endl;

            for ( int l=0; l<NUM_LAYER; ++l ) {
                for ( int v=0; v<NUM_VIEW; ++v ) {
                    std::cout << "layer " << l << " view " << v;
                    std::cout << " ToT " << ToT[0][l][v] << " " << ToT[1][l][v];
                    std::cout << " numHits: " << TkrNumHits[l][v];
                    if ( TkrNumHits[l][v] > 0 )
                        std::cout << " strips";
                    for ( int j=0; j<128 && TkrHits[l][v][j]>=0; ++j )
                        std::cout << " " << TkrHits[l][v][j];
                    std::cout << std::endl;
                }
            }

            std::cout << "TkrNumClus " << TkrNumClus << std::endl;
            std::cout << "TkrNumTracks " << TkrNumTracks << std::endl;
            if ( TkrNumTracks > 0 ) {
                std::cout << "Track 1 clusters:";
                int i = 0;
                for ( ; i<TkrNumClus && TkrTrk1Clusters[i]>=0; ++i )
                    std::cout << " " << TkrTrk1Clusters[i];
                std::cout << " (" << i << ")" << std::endl;
            }
        }

        for ( int l=0; l<NUM_LAYER; ++l )
            for ( int v=0; v<NUM_VIEW; ++v )
                for ( int j=0; j<TkrNumHits[l][v]; ++j )
                    TkrDigiHitMap[l][v]->Fill(TkrHits[l][v][j]);

        // the formula is: z = m * q + b, with q being x or y
        double Sq[2]  = { 0, 0 };
        double Sz[2]  = { 0, 0 };
        double Sqq[2] = { 0, 0 };
        double Szz[2] = { 0, 0 };
        double Sqz[2] = { 0, 0 };
        int n[2] = { 0, 0 };

        if ( TkrNumTracks > 0 ) {
            for ( int i=0; i<TkrNumClus && TkrTrk1Clusters[i]>=0; ++i ) {
                const int id = TkrTrk1Clusters[i];
                const int v    = TkrClusView[id];
                const double q = v ? TkrClusY[id] : TkrClusX[id];
                const double z = TkrClusZ[id];
                Sq[v]  += q;
                Sz[v]  += z;
                Sqq[v] += q * q;
                Szz[v] += z * z;
                Sqz[v] += q * z;
                ++n[v];
            }
            double m[2];
            double b[2];
            double r[2];
            // calculating the line parameters
            for ( int v=0; v<2; ++v ) {
                m[v] = (n[v]*Sqz[v]-Sq[v]*Sz[v]) / (n[v]*Sqq[v]-Sq[v]*Sq[v]);
                b[v] = ( Sz[v] - m[v] * Sq[v] ) / n[v];
                r[v] = ( n[v] * Sqz[v] - Sq[v] * Sz[v] )
                    / sqrt((n[v]*Sqq[v]-Sq[v]*Sq[v])*(n[v]*Szz[v]-Sz[v]*Sz[v]));
                if ( DEBUG )
                    std::cout << "view " << v << ":"
                              << " n " /* << std::setw(2) */ << n[v]
                              << " m " /* << std::setw(8) << std::fixed
                                          << std::setprecision(3) */ << m[v]
                              << " b " /* << std::setw(8) << std::fixed
                                          << std::setprecision(3) */ << b[v]
                              << " r " /* << std::setw(8) << std::fixed
                                          << std::setprecision(3) */ << r[v]
                              << std::endl;
            }
            // calculating the residuals
            double residual[128]; // Not all fields will be initialized!
            for ( int i=0; i<TkrNumClus && TkrTrk1Clusters[i]>=0; ++i ) {
                const int id = TkrTrk1Clusters[i];
                const int v    = TkrClusView[id];
                const double pos[3] = {TkrClusX[id],TkrClusY[id],TkrClusZ[id]};
                // q = ( z - b ) / m
                const double fitPos = ( pos[2] - b[v] ) / m[v];
                residual[id] = pos[v] - fitPos;
            }
            // do some printing
            if ( DEBUG ) {
                for ( int i=0; i<TkrNumClus && TkrTrk1Clusters[i]>=0; ++i ) {
                    const int id = TkrTrk1Clusters[i];
                    const int l    = TkrClusPlane[id];
                    const int v    = TkrClusView[id];
                    const double pos[3] = { TkrClusX[id], TkrClusY[id],
                                            TkrClusZ[id] };

                    std::cout << "cluster[" << std::setw(2) << id << "]:"
                        /* << std::setw(3) */ << l
                        /* << std::setw(2) */ << v;
                    for ( int j=0; j<3; ++j ) {
                        std::cout << " ";
                        if ( j == v )
                            std::cout << "(";
                        else
                            std::cout << " ";
                        std::cout /* << std::setw(8) << std::fixed
                                     << std::setprecision(3) */ << pos[j];
                        if ( j == v )
                            std::cout << ")";
                        else
                            std::cout << " ";
                    }
                    std::cout << " residual " /* << std::setw(6) */
                              << residual[id];
                    for ( int it=0; it<(int)(fabs(residual[id])/0.228); ++it )
                        std::cout << " <";
                    std::cout << std::endl;
                }
            }
        }
    }
 
    TH1D* TkrDigiHitMapProfile[NUM_LAYER][NUM_VIEW];

    for ( int l=0; l<NUM_LAYER; ++l ) {
        for ( int v=0; v<NUM_VIEW; ++v ) {
            TH1F* h = TkrDigiHitMap[l][v];

            int nBins = h->GetNbinsX();
            TString name = h->GetName();
            //                int mini = (int)h->GetMinimum();
            int maxi = (int)h->GetMaximum();

            TString s = TrayId(l, v);
            TkrDigiHitMapProfile[l][v] = new TH1D("TkrDigiHitMapProfile"+s,
                                             "profile of the TkrDigi hit map"+s,
                                                  maxi+1, -0.5, maxi+0.5);

            // the method to define dead/noisy strips has to be defined
            if ( 0 ) {
                if ( h->GetBinContent(0) != 0 )
                    std::cout << name << ": underflow bins" << std::endl;
                if ( h->GetBinContent(nBins+1) != 0 )
                    std::cout << name << ": overflow bins" << std::endl;
                if ( maxi <= 0 )
                    std::cout << name << ": is empty" << std::endl;
                else {
                    std::cout << name
                              << ": inefficient strips (<50% of maximum of "
                              << maxi << "):";
                    for ( int i=0; i<nBins; ++i) {
                        const int value = (int)h->GetBinContent(i+1);
                        if ( 2*value < maxi )
                            std::cout << " " << i;
                    }
                    std::cout << std::endl;
                }
            }

            for ( int i=0; i<nBins; ++i )
                TkrDigiHitMapProfile[l][v]->Fill((int)h->GetBinContent(i+1));
        }
    }

    myOutFile.Write();

}
