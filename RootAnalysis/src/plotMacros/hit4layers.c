// Author: Alicia Kavelaars (May 2001)
//
// This function goes creates 4  plots for 4 selected layers, 
// showing the hit map in each layer and creating a ps file.

void hit4layers(char* histfile) {

  // Read in the histogram file
  TFile *histFile = new TFile(histfile);
  TNtuple *ntuple1 = (TNtuple*)histFile->Get("TKRLayStrip");
  TNtuple *ntuple2 = (TNtuple*)histFile->Get("TKR_ELS");
  TNtuple *ntuple3 = (TNtuple*)histFile->Get("N3");
  TNtuple *ntuple4 = (TNtuple*)histFile->Get("N4");

  // Print mean and number of entries on histograms
  gStyle->SetOptStat(110);
  
  // Canvas 1
  TCanvas* c1= new TCanvas("c1","Canvas 1",800,600);

  // Create our pads
  TPad* pad11 = new TPad("pad11","", 0.01, 0.48, 0.49, 0.94);
  TPad* pad12 = new TPad("pad12","", 0.51, 0.48, 0.99, 0.94);
  TPad* pad13 = new TPad("pad13","", 0.01, 0.01, 0.49, 0.47);
  TPad* pad14 = new TPad("pad14","", 0.51, 0.01, 0.99, 0.47);
  pad11->Draw(); pad12->Draw(); pad13->Draw(); pad14->Draw();
  pad11->SetGrid(); pad12->SetGrid();
  pad13->SetGrid(); pad14->SetGrid();

  // Title for our pads
 
  title1 = new TPaveLabel(0.2,0.95,0.8,0.99,"Number of hits");
  title1->SetFillColor(0);
  title1->Draw();


  // Define our histograms

  //****************************USER*******************************
  //Set here the title with the correct layers chosen 
  //Last two parameters set the # strips map 

  TH1F *h0 = new TH1F("StripsLay0","Layer 25",30,0,2000);
  TH1F *h1 = new TH1F("StripsLay1","Layer 24",30,0,2000);
  TH1F *h2 = new TH1F("StripsLay2","Layer 1",30,0,2000);
  TH1F *h3 = new TH1F("StripsLay3","Layer 0",30,0,2000);
 
  // Set the axis Labels
  
  h0->SetXTitle("Strip ID");
  h0->SetYTitle("Number of events/bin");
  h1->SetXTitle("Strip ID");
  h1->SetYTitle("Number of events/bin");
  h2->SetXTitle("Strip ID");
  h2->SetYTitle("Number of events/bin");
  h3->SetXTitle("Strip ID");
  h3->SetYTitle("Number of events/bin");
 
  //****************************USER*******************************
  //Set here the correct # of layers chosen
  
  // pad 11
  pad11->cd();
  //  pad11->SetLogy();
  ntuple1->Draw("strip>>StripsLay0","layer==25");
  h0->SetFillColor(38);
  h0->GetXaxis()->CenterTitle();
  h0->GetYaxis()->CenterTitle();
  h0->Draw();

  //pad 12
  pad12->cd();
  //  pad12->SetLogy();
  ntuple1->Draw("strip>>StripsLay1","layer==24");
  h1->SetFillColor(38);
  h1->GetXaxis()->CenterTitle();
  h1->GetYaxis()->CenterTitle();
  h1->Draw();

  //pad 13
  pad13->cd();
  //  pad13->SetLogy();
  ntuple1->Draw("strip>>StripsLay2","layer==1");
  h2->SetFillColor(38);
  h2->GetXaxis()->CenterTitle();
  h2->GetYaxis()->CenterTitle();
  h2->Draw();

  //pad 14
  pad14->cd();
  //  pad14->SetLogy();
  ntuple1->Draw("strip>>StripsLay3","layer==0");
  h3->SetFillColor(38);
  h3->GetXaxis()->CenterTitle();
  h3->GetYaxis()->CenterTitle();
  h3->Draw();

  //*****************************USER*******************************
  //This is the path to be set for the ps file to be created at 

  c1->Print("e:/Root_Files/Run/strips.ps");
  // c1->update();
}






