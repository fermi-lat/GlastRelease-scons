// Author: Alicia Kavelaars (May 2001)
//
// This function goes through a loop to create plots for each layer and
// show the map of TOT, creating a ps file.


void totmap(char* histfile) {

  gROOT->Reset();
  gStyle->SetOptStat(1100);

  gStyle->SetTitleColor(0); // Title box background Color
  gStyle->SetTitleH(0.11);  // Franction of title box hight
  gStyle->SetTitleW(0.35);  // Franction of title box width
  gStyle->SetTitleX(0.1);   // X position of the title box from left
  gStyle->SetTitleY(1.00);  // Y position of the title box from bottom

  TCanvas *c1 = new TCanvas("c1","BFEM TOT MAP",0,0,700,850);

  c1->SetTitle("BFEM TOT MAP");
  c1->GetFrame()->SetFillColor(46);
  c1->SetFillColor(19);
  c1->SetBorderMode(28);

  TText *c1title = new TText(0.2,0.96,"BFEM TKR TOT MAP");

  c1title->SetTextSize(0.04);
  c1title->Draw();
 
  TFile *histFile = new TFile(histfile); 
  TNtuple *TKRLayHits = (TNtuple *)histFile->Get("TKRLayHitsTOT");
  TNtuple *TKRLayStrip = (TNtuple *)histFile->Get("TKRLayStrip");

  int maxpad = 26;
  int i;

  TPad *hpad[26];

  hpad[0]  = new TPad("pad0" ,"This is pad0" ,0.00,0.82,0.25,0.95,19);
  hpad[1]  = new TPad("pad1" ,"This is pad1" ,0.25,0.82,0.50,0.95,19);
  hpad[2]  = new TPad("pad2" ,"This is pad2" ,0.50,0.82,0.75,0.95,19);
  hpad[3]  = new TPad("pad3" ,"This is pad3" ,0.75,0.82,1.00,0.95,19);
  hpad[4]  = new TPad("pad4" ,"This is pad4" ,0.00,0.69,0.25,0.82,19);
  hpad[5]  = new TPad("pad5" ,"This is pad5" ,0.25,0.69,0.50,0.82,19);
  hpad[6]  = new TPad("pad6" ,"This is pad6" ,0.50,0.69,0.75,0.82,19);
  hpad[7]  = new TPad("pad7" ,"This is pad7" ,0.75,0.69,1.00,0.82,19);
  hpad[8]  = new TPad("pad8" ,"This is pad8" ,0.00,0.56,0.25,0.69,19);
  hpad[9]  = new TPad("pad9" ,"This is pad9" ,0.25,0.56,0.50,0.69,19);
  hpad[10] = new TPad("pad10","This is pad10",0.50,0.56,0.75,0.69,19);
  hpad[11] = new TPad("pad11","This is pad11",0.75,0.56,1.00,0.69,19);
  hpad[12] = new TPad("pad12","This is pad12",0.00,0.43,0.25,0.56,19);
  hpad[13] = new TPad("pad13","This is pad13",0.25,0.43,0.50,0.56,19);
  hpad[14] = new TPad("pad14","This is pad14",0.50,0.43,0.75,0.56,19);
  hpad[15] = new TPad("pad15","This is pad15",0.75,0.43,1.00,0.56,19);
  hpad[16] = new TPad("pad16","This is pad16",0.00,0.30,0.25,0.43,19);
  hpad[17] = new TPad("pad17","This is pad17",0.25,0.30,0.50,0.43,19);
  hpad[18] = new TPad("pad18","This is pad18",0.50,0.30,0.75,0.43,19);
  hpad[19] = new TPad("pad19","This is pad19",0.75,0.30,1.00,0.43,19);
  hpad[20] = new TPad("pad20","This is pad20",0.00,0.17,0.25,0.30,19);
  hpad[21] = new TPad("pad21","This is pad21",0.25,0.17,0.50,0.30,19);
  hpad[22] = new TPad("pad22","This is pad22",0.50,0.17,0.75,0.30,19);
  hpad[23] = new TPad("pad23","This is pad23",0.75,0.17,1.00,0.30,19);
  hpad[24] = new TPad("pad24","This is pad24",0.00,0.04,0.25,0.17,19);
  hpad[25] = new TPad("pad25","This is pad25",0.25,0.04,0.50,0.17,19);

  TH1F *nhits[26];
 
  //***********************USER************************************
  //Last two parameters of nhits[#] can be changed according to 
  //the # counts in each layer
    nhits[0]  = new TH1F("totb0", "Plane12-X(L25)",100,-0.5,35.5);
    nhits[1]  = new TH1F("tota1", "Plane12-Y(L24)",100,-0.5,35.5);
    nhits[2]  = new TH1F("totb2", "Plane11-Y(L23)",100,-0.5,35.5);
    nhits[3]  = new TH1F("tota3", "Plane11-X(L22)",100,-0.5,35.5);
    nhits[4]  = new TH1F("totb4", "Plane10-X(L21)",100,-0.5,35.5);
    nhits[5]  = new TH1F("tota5", "Plane10-Y(L20)",100,-0.5,35.5);
    nhits[6]  = new TH1F("totb6", "Plane9-Y(L19)",100,-0.5,35.5);
    nhits[7]  = new TH1F("tota7", "Plane9-X(L18)",100,-0.5,35.5);
    nhits[8]  = new TH1F("totb8", "Plane8-X(L17)",100,-0.5,35.5);
    nhits[9]  = new TH1F("tota9", "Plane8-Y(L16)",100,-0.5,35.5);
    nhits[10] = new TH1F("totb10","Plane7-Y(L15)",100,-0.5,35.5);
    nhits[11] = new TH1F("tota11","Plane7-X(L14)",100,-0.5,35.5);
    nhits[12] = new TH1F("totb12","Plane6-X(L13)",100,-0.5,35.5);
    nhits[13] = new TH1F("tota13","Plane6-Y(L12)",100,-0.5,35.5);
    nhits[14] = new TH1F("totb14","Plane5-Y(L11)",100,-0.5,35.5);
    nhits[15] = new TH1F("tota15","Plane5-X(L10)",100,-0.5,35.5);
    nhits[16] = new TH1F("totb16","Plane4-X(L9)",100,-0.5,35.5);
    nhits[17] = new TH1F("tota17","Plane4-Y(L8)",100,-0.5,35.5);
    nhits[18] = new TH1F("totb18","Plane3-Y(L7)",100,-0.5,35.5);
    nhits[19] = new TH1F("tota19","Plane3-X(L6)",100,-0.5,35.5);
    nhits[20] = new TH1F("totb20","Plane2-X(L5)",100,-0.5,35.5);
    nhits[21] = new TH1F("tota21","Plane2-Y(L4)",100,-0.5,35.5);
    nhits[22] = new TH1F("totb22","Plane1-Y(L3)",100,-0.5,35.5);
    nhits[23] = new TH1F("tota23","Plane1-X(L2)",100,-0.5,35.5);
    nhits[24] = new TH1F("totb24","Plane0-X(L1)",100,-0.5,35.5);
    nhits[25] = new TH1F("tota25","Plane0-Y(L0)",100,-0.5,35.5);
   


  hpad[0]->Draw();
  hpad[0]->cd();
  hpad[0]->SetLogy(1);
  hpad[0]->GetFrame()->SetFillColor(0);

  TCut l31cut = "layer==31";
  TCut l30cut = "layer==30";
  TCut l29cut = "layer==29";
  TCut l28cut = "layer==28";
  TCut l27cut = "layer==27";
  TCut l26cut = "layer==26";
  TCut l25cut = "layer==25";
  TCut l24cut = "layer==24";
  TCut l23cut = "layer==23";
  TCut l22cut = "layer==22";
  TCut l21cut = "layer==21";
  TCut l20cut = "layer==20";
  TCut l19cut = "layer==19";
  TCut l18cut = "layer==18";
  TCut l17cut = "layer==17";
  TCut l16cut = "layer==16";
  TCut l15cut = "layer==15";
  TCut l14cut = "layer==14";
  TCut l13cut = "layer==13";
  TCut l12cut = "layer==12";
  TCut l11cut = "layer==11";
  TCut l10cut = "layer==10";
  TCut l9cut = "layer==9";
  TCut l8cut = "layer==8";
  TCut l7cut = "layer==7";
  TCut l6cut = "layer==6";
  TCut l5cut = "layer==5";
  TCut l4cut = "layer==4";
  TCut l3cut = "layer==3";
  TCut l2cut = "layer==2";
  TCut l1cut = "layer==1";
  TCut l0cut = "layer==0";

 
    TKRLayHits->Draw("totb*.2>>totb0",l25cut);
    TKRLayHits->Draw("tota*.2>>tota1",l24cut);
    TKRLayHits->Draw("totb*.2>>totb2",l23cut);
    TKRLayHits->Draw("tota*.2>>tota3",l22cut);
    TKRLayHits->Draw("totb*.2>>totb4",l21cut);
    TKRLayHits->Draw("tota*.2>>tota5",l20cut);
    TKRLayHits->Draw("totb*.2>>totb6",l19cut);
    TKRLayHits->Draw("tota*.2>>tota7",l18cut);
    TKRLayHits->Draw("totb*.2>>totb8",l17cut);
    TKRLayHits->Draw("tota*.2>>tota9",l16cut);
    TKRLayHits->Draw("totb*.2>>totb10",l15cut);
    TKRLayHits->Draw("tota*.2>>tota11",l14cut);
    TKRLayHits->Draw("totb*.2>>totb12",l13cut);
    TKRLayHits->Draw("tota*.2>>tota13",l12cut);
    TKRLayHits->Draw("totb*.2>>totb14",l11cut);
    TKRLayHits->Draw("tota*.2>>tota15",l10cut);
    TKRLayHits->Draw("totb*.2>>totb16",l9cut);
    TKRLayHits->Draw("tota*.2>>tota17",l8cut);
    TKRLayHits->Draw("totb*.2>>totb18",l7cut);
    TKRLayHits->Draw("tota*.2>>tota19",l6cut);
    TKRLayHits->Draw("totb*.2>>totb20",l5cut);
    TKRLayHits->Draw("tota*.2>>tota21",l4cut);
    TKRLayHits->Draw("totb*.2>>totb22",l3cut);
    TKRLayHits->Draw("tota*.2>>tota23",l2cut);
    TKRLayHits->Draw("totb*.2>>totb24",l1cut);
    TKRLayHits->Draw("tota*.2>>tota25",l0cut);
 
 

  for(i=0;i<maxpad;i++) {
    hpad[i]->Draw();
    hpad[i]->cd();
    hpad[i]->SetLogy(0);
    hpad[i]->GetFrame()->SetFillColor(10);
    hpad[i]->SetGrid();   
    nhits[i]->SetFillColor(46);
    nhits[i]->SetLabelSize(0.07,"X");
    nhits[i]->SetLabelSize(0.07,"Y");
    //****************************USER**********************************
    //range is set on each layer, this is a default
    nhits[i]->GetXaxis()->SetRange(0,500);
    //****************************USER**********************************
    //SetMaximum sets the Y-axis max range upon # of events selected;
    // i.e., 10 if Go(10), 50 if Go(1000)
    nhits[i]->SetMaximum(10);
    nhits[i]->SetXTitle("TOT(us)");
    nhits[i]->SetYTitle("Nevents");
    nhits[i]->SetTitleOffset(1.3,"X");
    nhits[i]->SetTitleOffset(1.3,"Y");
    nhits[i]->SetTitleSize(0.09,"X");
    nhits[i]->SetTitleSize(0.09,"Y");
    nhits[i]->Draw();
    c1->cd();
  }
  //*****************************USER*******************************
  //This is the path to be set for the ps file to be created at 
  c1->Print("e:/Root_Files/Run/HitsTOT.ps");
}













