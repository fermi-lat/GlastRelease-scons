// Alicia Kavelaars (June 2001)
//
// Prints all histograms in a psfile in separate pages 

{
    
    //gROOT->Reset();
    
    TFile *Hist = new TFile("Histograms.root");
    
    TPostScript *ps1 = new TPostScript("reconhistograms.eps",111);
    // ps1.Range(16,24); // width,height in cm
    
    
    //gStyle->SetOptStat(111111);
    
    gStyle->SetTitleColor(0); // Title box background Color
    gStyle->SetTitleH(0.11);  // Fraction of title box height
    gStyle->SetTitleW(0.70);  // Fraction of title box width
    gStyle->SetTitleX(0.1);   // X position of the title box from left
    gStyle->SetTitleY(1.00);  // Y position of the title box from bottom
    gStyle->SetCanvasColor(19);
    gStyle->SetFrameFillColor(19);//21 is best
    gStyle->SetTitleFont(2.8);
    gStyle->SetPadColor(19);
    gStyle->SetOptTitle(1);
    gStyle->SetTitleXSize(1);
    gStyle->SetStatFont(2.9);
    gStyle->SetStatStyle(1111);
    //  gStyle->SetStatFormat(2);
    gStyle->SetStatTextColor(4);
    gStyle->SetStatColor(10);
    gStyle->SetStatH(0.20);  // Franction of title box hight
    gStyle->SetStatW(0.20);  // Franction of title box width
    gStyle->SetStatX(1.0);   // X position of the title box from left
    gStyle->SetStatY(0.9);  // Y position of the title box from bottom
    
    
    
    TCanvas *c1 = new TCanvas("c1","RH1",0,0,800,800);
    TCanvas *c2 = new TCanvas("c2","RH2",0,0,800,800);
    TCanvas *c3 = new TCanvas("c3","RH3",0,0,800,800);
    TCanvas *c4 = new TCanvas("c4","RH4",0,0,800,800);
    TCanvas *c5 = new TCanvas("c5","RH5",0,0,800,800);
    TCanvas *c6 = new TCanvas("c6","RH6",0,0,800,800);
    
    c1->SetTitle("ReconHistograms");
    c1->GetFrame()->SetFillColor(46);
    c1->SetFillColor(19);
    c1->SetBorderMode(19);
    
    
    c1.Divide(2,3);
    
    // Pick up the desired histograms by name
    c1 -> cd(1);
    TH1F *h1 = (TH1F*)gDirectory->Get("CENTERSTRIP");
    h1->Draw();
    c1 -> cd(2);
    TH1F *h2 = (TH1F*)gDirectory->Get("NCLUSLAY");
    h2->Draw();
    c1 -> cd(3);
    TH1F *h3 = (TH1F*)gDirectory->Get("NCLUSTRIP");
    h3->Draw();
    c1 -> cd(4);
    TH1F *h4 = (TH1F*)gDirectory->Get("CLUSPOS");
    h4->Draw();  
    c1 -> cd(5);
    TH1F *h5 = (TH1F*)gDirectory->Get("NCLUSXY");
    h5->Draw();
    c1 -> cd(6);
    TH1F *h6 = (TH1F*)gDirectory->Get("NCLUSZ");
    h6->Draw();
    
    c1.Update();
    ps1->NewPage();
    c2.Divide(2,3);
    
    // Pick up the next set of histograms by name
    c2 -> cd(1); 
    TH1F *h7 = (TH1F*)gDirectory->Get("CLUSID");
    h7->Draw();
    c2 -> cd(2); 
    TH1F *h8 = (TH1F*)gDirectory->Get("CHI2");
    h8->Draw();
    c2 -> cd(3); 
    TH1F *h9 = (TH1F*)gDirectory->Get("EDET");
    h9->Draw();
    c2 -> cd(4);
    TH1F *h10 = (TH1F*)gDirectory->Get("TRACKTYPE");
    h10->Draw();
    c2 -> cd(5); 
    TH1F *h11 = (TH1F*)gDirectory->Get("EINPUT");
    h11->Draw();
    c2 -> cd(6); 
    TH1F *h12 = (TH1F*)gDirectory->Get("NCLUSZ");
    h12->Draw();
    
    c2.Update();
    ps1->NewPage();
    c3.Divide(2,3);
    
    // Pick up the next set of histograms by name
    c3 -> cd(1);  
    TH1F *h13 = (TH1F*)gDirectory->Get("QUALITY");
    h13->Draw();
    c3 -> cd(2);  
    TH1F *h14 = (TH1F*)gDirectory->Get("TRACKTYPE");
    h14->Draw();
    c3 -> cd(3); 
    TH1F *h15 = (TH1F*)gDirectory->Get("RESIDUAL");
    h15->Draw();
    c3 -> cd(4);  
    TH1F *h16 = (TH1F*)gDirectory->Get("TRACKID");
    h16->Draw();
    c3 -> cd(5);  
    TH1F *h17 = (TH1F*)gDirectory->Get("FIRSTLAYER");
    h17->Draw();
    c3 -> cd(6);  
    TH1F *h18= (TH1F*)gDirectory->Get("NUMCLUSTERS");
    h18->Draw();
    
    c3.Update();
    ps1->NewPage();
    
    c4.Divide(2,3);
    // Pick up the next set of hisgtograms for the next page
    c4 -> cd(1);  
    TH1F *h19= (TH1F*)gDirectory->Get("HITCLUS");
    h19->Draw();
    c4 -> cd(2);  
    TH1F *h20= (TH1F*)gDirectory->Get("HITLOCATOR");
    h20->Draw();
    c4 -> cd(3);  
    TH1F *h21= (TH1F*)gDirectory->Get("HITRESIDUAL");
    h21->Draw();
    c4 -> cd(4); 
    TH1F *h22= (TH1F*)gDirectory->Get("HITCHI2");
    h22->Draw();
    c4 -> cd(5); 
    TH1F *h23= (TH1F*)gDirectory->Get("ELOGREC");
    h23->Draw();
    c4 -> cd(6); 
    TH1F *h24= (TH1F*)gDirectory->Get("ETOTREC");
    h24->Draw();
    
    
    c4.Update();
    ps1->NewPage();
    
    c5.Divide(2,3);
    c5 -> cd(1);  
    TH1F *h25= (TH1F*)gDirectory->Get("XSLOPE");
    h25->Draw();
    c5 -> cd(2);  
    TH1F *h26= (TH1F*)gDirectory->Get("CSIALPHA");
    h26->Draw();
    c5 -> cd(3); 
    TH1F *h27= (TH1F*)gDirectory->Get("CSILAMBDA");
    h27->Draw();
    c5 -> cd(4); 
    TH1F *h28= (TH1F*)gDirectory->Get("CSISTART");
    h28->Draw();
    c5 -> cd(5);  
    TH1F *h29= (TH1F*)gDirectory->Get("NLOGSREC");
    h29->Draw();
    c5 -> cd(6);  
    TH1F *h30= (TH1F*)gDirectory->Get("ECORR");
    h30->Draw();
    
    
    c5.Update();
    ps1->NewPage();
    
    c6.Divide(2,3);
    c6 -> cd(1);  
    TH1F *h33= (TH1F*)gDirectory->Get("ELAYER");
    h33->Draw();
    c6 -> cd(2);  
    TH1F *h34= (TH1F*)gDirectory->Get("EFIT");
    h34->Draw();
    TH1F *h35= (TH1F*)gDirectory->Get("ECLUSID");
    h35->Draw();
    c6 -> cd(3);  
    TH1F *h36= (TH1F*)gDirectory->Get("ECHISQ");
    h36->Draw();
    c6 -> cd(4); 
    TH1F *h37= (TH1F*)gDirectory->Get("ESUM");
    h37->Draw();
    c6 -> cd(5);
    TH1F *h38= (TH1F*)gDirectory->Get("ECLUSID");
    h38->Draw();
    
    c6.Update();
    
    ps1->Close();
    
    //*****************************USER***********************************
    //This is the path to be set for the ps file to be created at 
    //c1->Print("RH1.ps");
    //c2->Print("RH2.ps");
    //c3->Print("RH3.ps");
    //c4->Print("RH4.ps");
    //c5->Print("RH5.ps");
    //c6->Print("RH6.ps");
    // c7->Print("RH7.ps");
    
    
    
}



















