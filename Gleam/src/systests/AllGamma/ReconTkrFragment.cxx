//
// This file purports to fill some basic histograms for system tests
// Basically, it will provide some summary track and vertex information
// 
// Warranty expires upon first loading of macro into root.
//
// September 10, 2002  Tracy Usher
//
void RootTreeAnalysis::ReconTkrHistDefine() 
{
    
    histFile->cd();

    TH1F* dummy;

    dummy = new TH1F("TKRNUMTRACKS",   "Tkr number fit tracks",  20,    0.,   20.);
    dummy = new TH1F("TKRNUMVERTICES", "Tkr number vertices",    10,    0.,   10.);

    dummy = new TH1F("TKRTRKNHITS",    "Track Number of hits",   40,    0.,   40.);
    dummy = new TH1F("TKRTRKQUALITY",  "Track Quality",          80,    0.,   80.);
    dummy = new TH1F("TKRTRKSTARTZ",   "Track start Z",         100,   50.,  650.);
    dummy = new TH1F("TKRTRKSTARTX",   "Track start X",         100, -500.,  500.);
    dummy = new TH1F("TKRTRKSTARTY",   "Track start Y",         100, -500.,  500.);
    dummy = new TH1F("TKRTRKSLOPEX",   "Track slope X",         100,   -5.,    5.);
    dummy = new TH1F("TKRTRKSLOPEY",   "Track slope Y",         100,   -5.,    5.);
    dummy = new TH1F("TKRTRKENERGY",   "Track Energy",          200,    0., 2000.);

    dummy = new TH1F("TKRHITRADLEN",   "Track RadLen/plane",    100,    0.,    1.);
    dummy = new TH1F("TKRHITACTDIS",   "Active Distance",       100, -100.,  100.);
    dummy = new TH1F("TKRHITPLANEZ",   "Track Hit Plane Z",     600,   50.,  650.);

    dummy = new TH1F("TKRVTXQUALITY",  "Vertex Quality",        100,    0.,   0.5);
    dummy = new TH1F("TKRVTXSTARTZ",   "Vertex Z position",     100,   50.,  650.);
    dummy = new TH1F("TKRVTXNUMTRK",   "# Fit Tracks/Vertex",     5,    0,      5);

    return;
}


TkrTrack* findTrack(TObjArray* trackCol, int idx)
{
    TkrTrack* track     = 0;
    int       numTracks = trackCol->GetEntries();

    while(numTracks--)
    {
        TkrTrack* test = (TkrTrack*)trackCol->At(numTracks);

        if (test->getId() == idx) 
        {
            track = test;
            break;
        }
    }

    return track;
}


void RootTreeAnalysis::ReconTkr() 
{
    TkrRecon*   tkrEvt      = rec->getTkrRecon();
    TObjArray*  vtxCol      = tkrEvt->getVertexCol();
    TObjArray*  trkCol      = tkrEvt->getTrackCol();

    Int_t       numVertices = vtxCol->GetEntries();
    Int_t       numTracks   = trkCol->GetEntries();

    ((TH1F*)GetObjectPtr("TKRNUMVERTICES"))->Fill(numVertices);
    ((TH1F*)GetObjectPtr("TKRNUMTRACKS"))->Fill(numTracks);

    for (int idx = 0; idx < numTracks; idx++)
    {
        TkrTrack* track = (TkrTrack*)trkCol->At(idx);

        Int_t      numHits = track->getNumHits();
        Double_t   quality = track->getQuality();
        Double_t   trkZ    = track->getTrackParZ(TkrTrack::FIRSTHIT);
        Double_t   trkE    = track->getKalEnergy();
        const TkrParams* trkPrms = track->getTrackPar(TkrTrack::FIRSTHIT);

        ((TH1F*)GetObjectPtr("TKRTRKNHITS"))->Fill(numHits);
        ((TH1F*)GetObjectPtr("TKRTRKQUALITY"))->Fill(quality);
        ((TH1F*)GetObjectPtr("TKRTRKSTARTZ"))->Fill(trkZ);
        ((TH1F*)GetObjectPtr("TKRTRKSTARTX"))->Fill(trkPrms->getXPos());
        ((TH1F*)GetObjectPtr("TKRTRKSTARTY"))->Fill(trkPrms->getYPos());
        ((TH1F*)GetObjectPtr("TKRTRKSLOPEX"))->Fill(trkPrms->getXSlope());
        ((TH1F*)GetObjectPtr("TKRTRKSLOPEY"))->Fill(trkPrms->getYSlope());
        ((TH1F*)GetObjectPtr("TKRTRKENERGY"))->Fill(trkE);

	for(int planeIdx=0; planeIdx < track->getNumHits(); planeIdx++)
        {
            TkrHitPlane plane = track->getHitPlane(planeIdx);

            Double_t radLen = plane.getRadLen();
            Double_t planeZ = plane.getZplane();

            ((TH1F*)GetObjectPtr("TKRHITRADLEN"))->Fill(radLen);
            ((TH1F*)GetObjectPtr("TKRHITACTDIS"))->Fill(plane.getActiveDist());
            ((TH1F*)GetObjectPtr("TKRHITPLANEZ"))->Fill(planeZ);
        }
    }

    for (idx = 0; idx < numVertices; idx++)
    {
        TkrVertex* vertex = (TkrVertex*)vtxCol->At(idx);

        // Temporary to avoid bad vertices currently in TDS/PDS
        if (vertex->getEnergy() > 0. && vertex->getEnergy() < 500000.)
        {
            TVector3 vtxPos = vertex->getPosition();

            ((TH1F*)GetObjectPtr("TKRVTXQUALITY"))->Fill(vertex->getQuality());
            ((TH1F*)GetObjectPtr("TKRVTXSTARTZ"))->Fill(vtxPos.Z());
            ((TH1F*)GetObjectPtr("TKRVTXNUMTRK"))->Fill(vertex->getNumTracks());

/* Future expansion
            Double_t locQual = 0.;
            for(int jdx = 0; jdx < vertex->getNumTracks(); jdx++)
            {
                TkrTrack* track = findTrack(tkrEvt->getTrackCol(), vertex->getTrackId(jdx));

                locQual += track->getQuality();
            }

            locQual /= vertex->getNumTracks();

            if (locQual > vtxQual)
            {
                vtxQual   = locQual;
                bestIdx   = idx;
                foundBest = 1;
            }
*/
	}
    }
    return;
};
