// $Header$
// aeffmacro.C
//   plot aeff as function  of logE and cosT
//
// Usage: 
// > cd ${IRFANALYSIS}/data
// > root
// root[] .x ../macro/aeffmacro.C 
// ______________________________________________

{
gROOT->Reset();

TFile *f = new TFile("aeffTable.root");
f->ls();
TH2F* ha =(TH2F*)f->Get("lecta");
TH2F* hf =(TH2F*)f->Get("lectf");
TH2F* hb =(TH2F*)f->Get("lectb");

  TCanvas *ca = new TCanvas();
  ha->Draw("surf1");
  TCanvas *cf = new TCanvas();
  hf->Draw("surf1");
  TCanvas *cb = new TCanvas();
  hb->Draw("surf1");

Stat_t statf[7];
((TH2*)hf)->GetStats(statf);
 statf[0];


 TCanvas *c = new TCanvas();
 c->Divide(2,2);
 c->cd(2);
 ha->Draw("surf1");
 c->cd(3);
 hf->Draw("surf1");
 c->cd(4);
 hb->Draw("surf1");
 c->cd(1);
  }
