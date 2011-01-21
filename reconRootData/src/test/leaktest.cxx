// Test routine to make sure there are no memory leaks when creating
// and deleting reconRootData objects.
// To Run:
// 1) Make sure reconRootData.dll is in your ROOT library path
//    a) This can be done by either modifying your local .rootrc file 
//       Unix.*.Root.DynamicPath:    $(ROOTANALYSIS)/lib
//       WinNT.*.Root.DynamicPath:   $(ROOTANALYSIS)/lib
//   OR
//    b) Copy reconRootData.so (or .dll) into the directory from where you start ROOT.
// 2) You want to make sure that memory statistics are kept in ROOT by modifying your
//    local .rootrc file, and setting:
//    Root.MemStat:            1
//    Root.ObjectStat:         1
// 3) Start ROOT
// 4) At the ROOT prompt, type:  ".x leaktest.cxx" 

{
    UInt_t numEvents = 500;
    UInt_t numXtals = 100;
    UInt_t runNum = 1;
    const UInt_t numClusters = 20;
    const UInt_t numTracks = 15; 
    const UInt_t numVertices = 11;

    gObjectTable->Print();
    
    gSystem->Load("libcommonRootData.so");
    gSystem->Load("libreconRootData.so");

    TFile *file =  new TFile("recon.root", "RECREATE");
    TTree *t = new TTree("Recon", "Recon");
    ReconEvent *ev = new ReconEvent();
    t->Branch("ReconEvent", "ReconEvent", &ev, 64000, 1);
    
    gObjectTable->Print();
    //TkrCandTrack *candTrack;
    TkrTrack *track;
    TkrVertex *vertex;


    Int_t ievent, ixtal;
    for (ievent = 0; ievent < numEvents; ievent++) {
        Float_t randNum = gRandom->Rndm();
        Float_t f = Float_t (ievent);

        // Create AcdRecon object
        AcdRecon *acdRec = new AcdRecon();
        Double_t energy = f;
        Int_t count = 5;
        Double_t gDoca = f*randNum;
        Double_t doca = randNum;
        Double_t actDist = 2.*f;
        AcdId minDocaId(0, 0, 3, 2);
        acdRec->addRowDoca(randNum);
        acdRec->addRowDoca(f);
        acdRec->addRowActDist(randNum);
        acdRec->addRowActDist(f);
                acdRec->addEnergy(f*2.);
                acdRec->addEnergy(f);
                acdRec->addId(AcdId(0, 0, 1, 2));
                acdRec->addId(AcdId(0, 1, 3, 4));
        acdRec->initialize(energy, count, gDoca, doca, actDist, minDocaId);

                // Create CalRecon object
        CalRecon *calRec = new CalRecon();
        calRec->initialize();

        UInt_t icluster;
        for (icluster = 0; icluster < numClusters; icluster++ ) {
            CalCluster *cluster = new CalCluster();
            cluster->addEnergy(15.0);
            cluster->addEnergy(22.0);
            cluster->addPos(TVector3(f, f, f));
            cluster->addPos(TVector3(randNum*f, randNum*2.*f, randNum*3.*f));
            cluster->addRms(TVector3(randNum, randNum*5., randNum*10.));
            cluster->addRms(TVector3(f+1., f+2., f+3.));
            cluster->addRms(TVector3(14., 22., 44.));
            TVector3 calDir(randNum, randNum*3, randNum*5);
            Double_t eLeak = randNum;
            Double_t rmsLong = randNum*f;
            Double_t rmsTrans = randNum*2.*f;
            Double_t transOffset = randNum*3.*f;
            Double_t eCorrect = randNum*4.*f;
            cluster->initialize(eLeak, eCorrect, rmsLong, rmsTrans, calDir, transOffset);
            
            Double_t fitEnergy = f;
            Double_t chi2 = f*randNum;
            Double_t fStart = randNum*randNum;
            Double_t fitAlpha = f*f;
            Double_t fitLambda = f+f;
            cluster->initProfile(fitEnergy, chi2, fStart, fitAlpha, fitLambda);
            
            calRec->addCalCluster(cluster);
        }
        
        for (ixtal = 0; ixtal < numXtals; ixtal ++) {
            CalXtalRecData* xtal = new CalXtalRecData();
            CalXtalId id;
            id.init(1, 2, 3);
            xtal->initialize(CalXtalId::BESTRANGE, id);
            CalRangeRecData rec(CalXtalId::LEX8, randNum*f, CalXtalId::HEX8, randNum*4.0);
            TVector3 pos(4.5, 7.5, 8.5);
            rec.initialize(pos);
            xtal->addRangeRecData(rec);
            calRec->addXtalRecData(xtal);
        }

        // Create TkrRecon object
        TkrRecon *tkrRec = new TkrRecon();
        tkrRec->initialize();

        for (icluster = 0; icluster < numClusters; icluster++ ) 
        {
            UInt_t towerX = icluster % 4;
            UInt_t towerY = (icluster / 4) % 4;
            UInt_t tray   = icluster % 15;
            commonRootData::TkrId tkrId(towerX, towerY, tray, true, commonRootData::TkrId::eMeasureX);
            UInt_t strip0 = 4;
            UInt_t stripf = 10;
            TVector3 pos(f, 2.*f, 3.*f);
            UInt_t rawTot = 50;
            Double_t tot =  f;
            UInt_t flag = 1;
            UInt_t tower = 8;
            UInt_t nBad  = 0;
            TkrCluster *tcluster = new TkrCluster(tkrId, strip0, stripf, pos, rawTot, tot, flag, nBad);
            tkrRec->addCluster(tcluster);
        }

        TkrTrack *track;
        for (int itrack=0; itrack < numTracks; itrack++) 
        {
            track = new TkrTrack();
            UInt_t xgaps = itrack * itrack;
            UInt_t ygaps = itrack+itrack;
            UInt_t x1st = ievent+itrack;
            UInt_t y1st = ievent+ievent;
            track->setNumXGaps(xgaps);
            track->setNumYGaps(ygaps);
            track->setNumXFirstGaps(x1st);
            track->setNumYFirstGaps(y1st);
            Double_t chiSq = f;
            Double_t chiSqSmooth = f*randNum;
            Double_t rms = randNum*randNum;
            Double_t qual = f*f;
            Double_t e = f*2.0;
            Double_t ms = f*6.0;
            track->setChiSquareFilter(chiSq);
            track->setChiSquareSmooth(chiSqSmooth);
            track->setScatter(rms);
            track->setQuality(qual);
            track->setKalEnergy(e);
            track->setKalThetaMS(ms);
            tkrRec->addTrack(track);
        }

        UInt_t ivertex;
        TkrVertex *vertex;
        for (ivertex=0; ivertex < numVertices; ivertex++) {
            vertex = new TkrVertex();
            UInt_t towerX = ivertex % 4;
            UInt_t towerY = ivertex / 4;
            UInt_t tray   = ivertex % 15;
            commonRootData::TkrId tkrId(towerX, towerY, tray, true, commonRootData::TkrId::eMeasureX);
            vertex->setTkrID(tkrId);
            Double_t qual = f;
            vertex->setQuality(qual);
            Double_t energy = f*randNum;
            vertex->setEnergy(energy);

            TkrTrackParams vtxPar;
            vtxPar(1) = f;
            vtxPar(2) = f*f;
            vtxPar(3) = f*randNum;
            vtxPar(4) = randNum*randNum;

            vtxPar(1,1) = f;
            vtxPar(1,2) = f*2.;
            vtxPar(1,3) = f*3.;
            vtxPar(1,4) = f*4.;
            vtxPar(2,2) = f*5.;
            vtxPar(2,3) = f*6.;
            vtxPar(2,4) = f*7.;
            vtxPar(3,3) = f*8.;
            vtxPar(3,4) = f*9.;
            vtxPar(4,4) = f*10.;
            vertex->setParams(vtxPar);

            TVector3 pos(f, f*randNum, f*randNum*2.);
            vertex->setPosition(pos);
            TVector3 dir(randNum*2., randNum*3., randNum*4.);
            vertex->setDirection(dir);
            //vertex->addTrackId(ivertex+1);
            //vertex->addTrackId(ivertex+2);
            tkrRec->addVertex(vertex);
        }



        ev->initialize(ievent, runNum, tkrRec, calRec, acdRec);
        t->Fill();
        ev->Clear();
    }
    
    delete ev;
    
    printf("Here is the object table after creating events, storing them to file and deleting the objects\n");
    gObjectTable->Print();
    
    file->Write();
    file->Close();
    delete file;
    
}


