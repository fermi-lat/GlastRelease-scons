{
//   example of macro to read data from an ascii file and
//   create a root file with an histogram.
 
gROOT->Reset();

// This is the output root file
TFile *hfile = (TFile*)gROOT->FindObject("../../../st.root"); 
if (hfile) hfile->Close();
hfile = new TFile("../../../st.root","RECREATE","ST Root File");

// This is the input file produced by the ST test; note the path
FILE *fp1 = fopen("../../../test.dat","r");  

// This ntuple contains the relevant data
ntuple = new TNtuple("ntuple","Demo ntuple","Edep:Eloss:angle");

// The physical quantities
Float_t Edep, Eloss, angle;

Int_t ncols;
Int_t nlines = 0;


// We fill the ntuple
while (1) {
  ncols = fscanf(fp1,"%f %f %f",&Edep, &Eloss, &angle);
  if (ncols < 0) break;    
  ntuple->Fill(Edep, Eloss, angle);
  nlines++;
}

// We build some Histograms
ntuple->Draw("Edep>>hEdep");
ntuple->Draw("Eloss>>hEloss");
ntuple->Draw("angle>>hangle");

// We fit energy deposited as a Landau
hEdep->Fit("landau");

// We write everything
hfile->Write();

fclose(fp1);
}

