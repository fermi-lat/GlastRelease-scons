#include <string>
#include "tree_class.h"

static inline double sqr(double x){return x*x;}

int main(){

    using namespace std;
    /******************************************************
    Procedure:

    1. The ntuple file must be "fulltup.root"
    2. .L tree_class.C
    3. .x newbranches.C

    (If you want to use a differente file name you have to change the
    code here and in the definition of the class).

    ********************************************************/

    /**************************************************************
    Branches to be added to the ntuple:

    Tag_Id =  1  for Thin-Tkr  events
    2  for Thin-Vtx  events
    3  for Thick-Tkr events
    4  for Thick-Vtx events

    0  for events with bad energy. (IMgoodCalProb<0.5).


    Dir_Err = Angle (in degrees) between the Mc direction and the
    best reconstructed direction (Vtx or Trk1)

    ReconZDir = Best ZDir (Vtx or Trk1)

    ReconPhiDir = Best PhiDir (Vtx or Trk1)

    ****************************************************************/

//    gROOT->Reset();


    std::string path = ::getenv("file_root");
    TFile *f = new TFile( (path+"/fulltup.root").c_str(), "update");
    cout << " newbranches : read file "<< f->GetName() << endl;

    TTree *tree    = (TTree*)f->Get("1");
    Int_t nentries = (Int_t)tree->GetEntries();
    cout << "    get tree "<< tree->GetName() <<
            ", number entries "<<nentries<<endl;

    tree_class t(tree);
#if 0
    Float_t Tag_Id, DirErr, ReconZDir, ReconPhiDir;
    TBranch *tag = tree->Branch("Tag_Id", &Tag_Id, "Tag_Id/F");
    TBranch *err = tree->Branch("BestDirErr", &DirErr, "BestDirErr/F");
    TBranch *recon = tree->Branch("ReconZDir", &ReconZDir, "ReconZDir/F");
    TBranch *phi = tree->Branch("ReconPhiDir", &ReconPhiDir, "ReconPhiDir/F");
#endif
    float psf_scale_factor;
    TBranch *psfScaleFactor = tree->Branch("PSFscaleFactor", &psf_scale_factor, "PSFscaleFactor/F");

    for(int  i=0; i<nentries; i++){
        tree->GetEntry(i);
        t.GetEntry(i);
        //cout<<t.IMgoodCalProb;
#if 0
        if(t.IMgoodCalProb<0.5){
            Tag_Id      = 0.0;
            DirErr      = 0.0;
            ReconZDir   = 0.0;
            ReconPhiDir = 0.0;
        }
        else if (t.IMvertexProb<0.5||t.VtxAngle==0.0){
            if (t.Tkr1FirstLayer<12.0)
                Tag_Id=1.0;
            else
                Tag_Id=3.0;
            DirErr      = t.McTkr1DirErr;
            ReconZDir   = t.Tkr1ZDir;
            ReconPhiDir = t.Tkr1Phi;
        }
        else{
            if(t.Tkr1FirstLayer<12.0)
                Tag_Id=2.0;
            else
                Tag_Id=4.0;
            DirErr      = t.McDirErr;
            ReconZDir   = t.VtxZDir;
            ReconPhiDir = t.VtxPhi;
        }
#endif
        // scale factor for thin and thick Tkr layers
        psf_scale_factor=sqrt(sqr(t.Tkr1ThetaErr)+sqr(t.Tkr1PhiErr));
        if (t.Tkr1FirstLayer<12.0)
           psf_scale_factor   *= 2.5; 
        else psf_scale_factor *= 3.5;
#if 0
        tag->Fill();
        err->Fill();
        recon->Fill();
        phi->Fill();
#endif

        psfScaleFactor->Fill();
        //cout<<"event"<<i<<endl;
    }

    tree->Write("", TObject::kOverwrite);
    f->Close();
}

