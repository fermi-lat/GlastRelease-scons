{
    gSystem->Load("libdigiRootData.so");//{{{1
    gSystem->Load("libmcRootData.so");
    gSystem->Load("libreconRootData.so");
    gROOT->Reset();    
    gROOT->LoadMacro("/home/pol/Glast/Analyse/RootTreeAnalysis/Current/MuonAnalysis.cxx");//}}}1
    //gROOT->LoadMacro("/home/pol/Glast/Analyse/RootTreeAnalysis/Current/HistBuildtuple.cxx");

	//chains//{{{1
    char *digilist[]={//{{{2
        //"sym_pdl_digi2.root",
        "/disk1/data/calibData/muons_digi/ebf031004040459_digi.root",
        "/disk1/data/calibData/muons_digi/ebf031004060929_digi.root"
    };//}}}2
    
    //char *mclist[]={"", "","",""};
    char *reconlist[]={//{{{2
        //"sym_pdl_recon.root",
        "/disk1/data/calibData/muons_recon/ebf031004040459_recon.root",
        "/disk1/data/calibData/muons_recon/ebf031004060929_recon.root"
    };//}}}2

    int numFiles = sizeof(digilist)/sizeof(char*);//{{{2
    cout<<"Found "<<numFiles<<" file(s)"<<endl;

    char *digitreePath = "Digi";
    char *mctreePath = "Mc";
    char *recontreePath = "Recon";

    TChain *digichainedTree = new TChain(digitreePath);
    TChain *mcchainedTree;
    //TChain *reconchainedTree;
    TChain *reconchainedTree = new TChain(recontreePath);//}}}2
    for(int i=0; i < numFiles; i++){//{{{2
        printf("Adding file %s\n",digilist[i]);        
        digichainedTree->Add(digilist[i]);
        //printf("Adding file %s\n",mclist[i]);        
        //mcchainedTree->Add(mclist[i]);
        printf("Adding file %s\n",reconlist[i]);        
        reconchainedTree->Add(reconlist[i]);
    }//}}}2
    //}}}1
    short int job_id= 5;

    RootTreeAnalysis *m = new RootTreeAnalysis(//{{{1
					  digichainedTree,
					  reconchainedTree,
					  mcchainedTree, 
            "try.root",
            job_id);//}}}1

    //m.Rewind();
    //m.setjobid(job_id);//Pdl
    //m.Go(20000);
    //m.FitPdl();

    //m.Rewind();
    //m.setjobid(1); //Correlated Pdl
    //m.Go(120000);
    //m.FitCorrPdl();

    //m.Rewind();
    //m.setjobid(2); //Correlated Pdl, cut on muon path
    //m.Go(120000);
    //m.FitCorrPdl("Clean"); 

    //m.Rewind();
    //m.setjobid(3); //Log(POS/NEG)
    //m.Go(120000);
    //m.FitLogPoverN();

    //m.setjobid(4); //Gain
    //m.Go(120000);
    //m.FitGain(10);

    m.Rewind();
    m.setjobid(5); //Muon signal vs position
    m.Go(120000);
    //m.FitTapering();
}
