#ifndef TreeMaker_cxx
#define TreeMaker_cxx 1

#include "TreeMaker.h"
#include <vector>
#include <bitset>
#include <iomanip>
#include <map>
#include <utility>

UInt_t digiEventId, reconEventId, mcEventId;
UInt_t digiRunNum, reconRunNum, mcRunNum;

const TString TS = "";
const int MaxNumLayers=36;
const Int_t NumGTCC = 8;
const Int_t NumGTRC = 9;

/*
 * The order of read-out of the GTCC's (LAT-TD-00605 p.106
 * http://www-glast.slac.stanford.edu/IntegrationTest/ONLINE/docs/TEM.pdf):
 * 6, 3, 7, 2, 5, 0, 4, 1
 */
const int positionOfGTCC[8] /* in the TCloneArray */ = { 5, 7, 3, 1, 6, 4, 0, 2 };

#define DEBUG 0
#define DIAGN 1

// returns the position of GTCC i in the TCloneArray m_TkrDiagnostics of the digi file
int indexToGTCC ( int i ) { return positionOfGTCC[i]; }

// returns the index of a GTCC/GTRC pair in the TkrDiagnostics array
int getIndex ( int GTCC, int GTRC ) { return GTCC*NumGTRC+GTRC; }

// returns the index of a layer/ROC pair in the TkrDiagnostics array. Left is 0, right is 1.
int getIndex ( const TString layerName, const bool side ) {
    // assume e.g. "LayerX12"
    static const int first = 5;
    char view = layerName[first];
    int num = atoi(layerName(first+1,layerName.Length()-first-1).Data());
    int GTRC = num / 2;
    int GTCC = 42;
    if ( view == 'X' )
        if ( num % 2 )            // odd
            GTCC = 6 + side;
        else                      // even
            GTCC = 3 - side;
    else                      // 'Y'
        if ( num % 2 )            // odd
            GTCC = 5 - side;
        else                      // even
            GTCC = 0 + side;
    
    //    std::cout << std::setw(8) << layerName << ' ' << view << num << " side=" << side << ' ' << GTCC << ' ' << GTRC << std::endl;

    return getIndex(GTCC, GTRC);
}

void TreeMaker::McData() 
{
}

/*
void TreeMaker::CreateDigiTree(Int_t numEvents)
{
}
*/

void TreeMaker::CreateTree(Int_t numEvents)
{    
  Int_t TrueToT0[MaxNumLayers];
  Int_t TrueToT1[MaxNumLayers];
  Int_t TrueTkrNumHits[MaxNumLayers];
  Int_t TrueTkrHits[MaxNumLayers][128];
  
  Int_t ToT0;//[MaxNumLayers];
  Int_t ToT1;//[MaxNumLayers];
  Bool_t TriggerReq0;
  Bool_t TriggerReq1;
  Int_t TkrNumHits;
  Int_t TkrHits[128];
  Int_t EventId;
  Int_t RunId;
  Int_t LevelOneTrigger;
  Double_t EbfTime;
  Int_t TkrDiagnostics[NumGTCC*NumGTRC]; // this should be sorted into the "layers", but right now we don't know the order
  Int_t TkrTotalNumHits;  
  ////////// clusters
  int TkrNumClus;
  int TkrClusPlane[128]; // !!
  int TkrClusView[128];
  double TkrClusX[128];
  double TkrClusY[128];
  double TkrClusZ[128];
  //////////  tracks
  
  int TkrNumTracks;
  int TkrTrk1Clusters[128];
  double TkrTheta[128];
  
  ////////// vertex //////////    
  
  int TkrNumVtx;
  double TkrVtxX, TkrVtxY, TkrVtxZ;
  
  //////////////////// Recon //////////////////////////////
  double Theta,Phi,ThetaXZ,ThetaYZ;
  
  TreeCollection = new TObjArray();
  
  // mc branches:
  if (mcTree) {
    mcTree->SetBranchStatus("*", 0);    // disable all branches
    // Activate desired branches...
    //        mcTree->SetBranchStatus("m_eventId", 1);
    //        mcTree->SetBranchStatus("m_particleCol", 1);
    //        mcTree->SetBranchStatus("m_runId", 1);        
    //        mcTree->SetBranchStatus("m_integratingHitCol", 1);        
    //        mcTree->SetBranchStatus("m_positionHitCol", 1);        
  }
  
  if (digiTree) {
    digiTree->SetBranchStatus("*", 0);  // disable all branches
    // activate desired brances
    digiTree->SetBranchStatus("m_tkr*", 1);
    digiTree->SetBranchStatus("m_eventId", 1); 
    digiTree->SetBranchStatus("m_runId", 1);
    if ( DIAGN ) {
	digiTree->SetBranchStatus("m_levelOneTrigger", 1);
	digiTree->SetBranchStatus("m_ebfTime*", 1);
	digiTree->SetBranchStatus("m_numTkrDiagnostics", 1);
    }
    for(int i=0;i<MaxNumLayers;i++)
      {
	TString TS="Layer";
	if(i % 2 == 0) 
	  TS+="X";
	else 
	  TS+="Y";
	TS+=i/2; 
	if ( DEBUG ) std::cout<<TS<<std::endl;
	TTree *Layer = new TTree(TS,TS);
	Layer->Branch("ToT0",&ToT0,"ToT0/I");
	Layer->Branch("ToT1",&ToT1,"ToT1/I");
	if ( DIAGN ) {
	    Layer->Branch("TriggerReq0",&TriggerReq0,"TriggerReq0/B");
	    Layer->Branch("TriggerReq1",&TriggerReq1,"TriggerReq1/B");
        }
	Layer->Branch("TkrNumHits",&TkrNumHits,"TkrNumHits/I");
	Layer->Branch("TkrHits",TkrHits,"TkrHits[TkrNumHits]/I");
	TreeCollection->Add((TObject*) Layer);
      }
    TTree *Header = new TTree("Header","Header");
    Header->Branch("EventId",&EventId,"EventId/I");
    Header->Branch("TkrTotalNumHits",&TkrTotalNumHits,"TkrTotalNumHits/I");
    Header->Branch("RunId",&RunId,"RunId/I");
    if ( DIAGN ) {
	Header->Branch("LevelOneTrigger",&LevelOneTrigger,"LevelOneTrigger/I");
	TString tag("TkrDiagnostics[");
	tag += NumGTCC * NumGTRC;
	tag += "]/I";
	Header->Branch("TkrDiagnostics", TkrDiagnostics, tag);
	Header->Branch("EbfTime",&EbfTime,"EbfTime/D");
    }
    TreeCollection->Add(Header);
    
  }
  if ( DEBUG ) std::cout<<TreeCollection->GetEntries()<<std::endl;
  
  if (reconTree) {
    reconTree->SetBranchStatus("*",0);  // disable all branches
    // activate desired branches
    reconTree->SetBranchStatus("m_cal*", 1);  
    reconTree->SetBranchStatus("m_tkr*", 1);
    reconTree->SetBranchStatus("m_acd*", 1);
    reconTree->SetBranchStatus("m_eventId", 1); 
    reconTree->SetBranchStatus("m_runId", 1);
    
    TTree *ReconTree = new TTree("Recon","Recon");
    ////////// Clusters: ////////////////////////////////////////
    ReconTree->Branch("TkrNumClus",&TkrNumClus,"TkrNumClus/I");
    ReconTree->Branch("TkrClusX",TkrClusX,"TkrClusX[TkrNumClus]/D");
    ReconTree->Branch("TkrClusY",TkrClusY,"TkrClusY[TkrNumClus]/D");
    ReconTree->Branch("TkrClusZ",TkrClusZ,"TkrClusZ[TkrNumClus]/D");
    ReconTree->Branch("TkrClusPlane",TkrClusPlane,"TkrClusPlane[TkrNumClus]/I");
    ReconTree->Branch("TkrClusView",TkrClusView,"TkrClusView[TkrNumClus]/I");
    /////////// Tracks: ////////////////////////////////////////
    ReconTree->Branch("TkrNumTracks",&TkrNumTracks,"TkrNumTracks/I");
    ReconTree->Branch("TkrTrk1Clusters",TkrTrk1Clusters,"TkrTrk1Clusters[TkrNumClus]/I");
    ReconTree->Branch("TkrTheta",TkrTheta,"TkrTheta[TkrNumTracks]/D");
    
    /////////// Vertex: ////////////////////////////////////////
    ReconTree->Branch("TkrNumVtx",&TkrNumVtx,"TkrNumVtx/I");
    ReconTree->Branch("TkrVtxX",&TkrVtxX,"TkrVtxX/D");
    ReconTree->Branch("TkrVtxY",&TkrVtxY,"TkrVtxY/D");
    ReconTree->Branch("TkrVtxZ",&TkrVtxZ,"TkrVtxZ/D");
    
    //////////Recon //////////
    ReconTree->Branch("Theta",&Theta,"Theta/D");
    ReconTree->Branch("Phi",&Phi,"Phi/D");
    ReconTree->Branch("ThetaXZ",&ThetaXZ,"ThetaXZ/D");
    ReconTree->Branch("ThetaYZ",&ThetaYZ,"ThetaYZ/D");
    TreeCollection->Add(ReconTree);

  }
  
  // determine how many events to process
  Int_t nentries = GetEntries();
  std::cout << "\nNum Events in File is: " << nentries << std::endl;
  Int_t curI = 0;
  Int_t nMax = TMath::Min(numEvents+m_StartEvent,nentries);
  if (m_StartEvent == nentries) {
    std::cout << " all events in file read" << std::endl;
    return;
  }
  if (nentries <= 0) return;
  
  // Keep track of how many bytes we have read in from the data files
  Int_t nbytes = 0, nb = 0;
  
  int target_fraction = 0;
  //////////////////////////////////////////////////
  // BEGINNING OF EVENT LOOP
  for (Int_t ievent=m_StartEvent; ievent<nMax; ievent++, curI=ievent) 
    {
      if ( DEBUG ) std::cout << "event loop counter: " << ievent << std::endl;
      int nTkrDigi = 0;
      // Reset some arrays and variables
      TkrTotalNumHits=0;
      
      for(int i=0;i<MaxNumLayers;i++) 
	{
	  TrueTkrNumHits[i]=0;
	  TrueToT0[i]=-1;
	  TrueToT1[i]=-1;
	}
      
      for ( int h=0; h<128; ++h )
	TkrTrk1Clusters[h] = -1;
      
      int fraction = (int)( (nMax-m_StartEvent)>0 ? 100.0*(ievent-m_StartEvent)/(nMax-m_StartEvent) : 0 );
      if ( fraction > target_fraction ) {
	target_fraction = fraction;
	if ( fraction >= 10 )
	  target_fraction = target_fraction/10*10 + 9;
	std::cout << fraction << "% complete: event "<< ievent << std::endl;
      }
      
      if (mc)
	mc->Clear();
      if (evt)
	evt->Clear();
      if (rec)
	rec->Clear();
      
      digiEventId = 0; reconEventId = 0; mcEventId = 0;
      digiRunNum = 0; reconRunNum = 0; mcRunNum = 0;
      
      nb = GetEvent(ievent);
      
      nbytes += nb;
      //////////////////////////////////////////////////
      // Monte Carlo ONLY analysis
      if (mc) {  // if we have mc data process it
          mcEventId = mc->getEventId();
          mcRunNum  = mc->getRunId();
          if ( DEBUG ) std::cout << "MC: run number: " << mcRunNum << "  event id: " << mcEventId << std::endl;
          //      McData();
      }
      
      // Digi ONLY analysis
      if (evt) {
	  digiEventId = evt->getEventId(); 
	  digiRunNum = evt->getRunId();
          if ( DEBUG ) std::cout << "DIGI: run number: " << digiRunNum << "  event id: " << digiEventId <<std::endl;
	  EventId = ievent+1;
	  RunId   = digiRunNum;
          if ( DEBUG ) std::cout << "run number: " << RunId << "  event id: " << EventId << std::endl;

	  if ( DIAGN ) { // we shouldn't go here if the digi is generated from mc
	      L1T l1t = evt->getL1T();
	      LevelOneTrigger = l1t.getTriggerWord() & 0x1f; // there are more bits set, but why?
              bool tkrTrigger = l1t.getTkr3InARow(); // indicates that 3-in-a-row occurred in the TKR
              // ***************************************************************
              // it is not clear if the files taken at Pisa in August and
              // September 2004 contain a correct L1T.  E.g., 040826180009.ldf
              // should have been taken with internal trigger (00100), but the
              // trigger word is 10010.
              // ***************************************************************
              // print some info about the trigger
              static bool tkrTriggerWarning = true;
              if ( DEBUG || tkrTriggerWarning ) { // print at least once
                  if ( !tkrTrigger )
                      std::cerr << "Trigger primitive:  External trigger!" << std::endl;
                  std::cout << "LevelOneTrigger: " << LevelOneTrigger << std::endl;
              }
              tkrTriggerWarning = false;

              // Ebf time
	      static Double_t EbfTimeStart = evt->getEbfTimeSec();
	      EbfTime = ( evt->getEbfTimeSec() - EbfTimeStart ) + evt->getEbfTimeNanoSec() * 1E-9;
	      if ( DEBUG ) std::cout << "EbfTime: " << EbfTime << std::endl;

	      // and now, for the TkrDiagnostics
              int numGTCCentries = evt->getTkrDiagnosticCol()->GetEntries();
              if ( NumGTCC != numGTCCentries ) {
                  static bool NumGTCCwarning = true;
		  if ( DEBUG || NumGTCCwarning ) // print at least once
                      std::cerr << "NumGTCC=" << NumGTCC
                                << " != evt->getTkrDiagnosticCol()->GetEntries()=" << numGTCCentries << std::endl;
                  NumGTCCwarning = false;
              }
	      for ( int GTCC=0; GTCC<NumGTCC; ++GTCC ) { 
                  std::bitset<NumGTRC> word = 0;
                  if ( tkrTrigger )
                      word = evt->getTkrDiagnostic(indexToGTCC(GTCC))->getDataWord();
                  if ( DEBUG ) std::cout << "TkrDiagnosticData[" << GTCC << "] " << word << ' ';
		  for ( int GTRC=NumGTRC-1; GTRC>=0; --GTRC ) {
                      TkrDiagnostics[getIndex(GTCC,GTRC)] = word[GTRC];
                      if ( DEBUG ) std::cout << ( TkrDiagnostics[getIndex(GTCC,GTRC)] > 0 ) ? 1 : 0;
		  }
		  if ( DEBUG ) std::cout << std::endl;
              }
          }
     	  //////////////////////////////////////////////////
	  // digitkr:
	  const TObjArray* tkrDigiCol = 0;
	  tkrDigiCol = evt->getTkrDigiCol();
	  //////////////////////////////////////////////////
	  if (tkrDigiCol)
	    {
	      nTkrDigi = tkrDigiCol->GetEntries();
	      if ( DEBUG ) std::cout << "num tkrDigi " << nTkrDigi << std::endl;
	      
	      TIter tkrDigiIter(tkrDigiCol);
	      TkrDigi* pTkrDigi = 0;
	      
	      
	      while ( (pTkrDigi=(TkrDigi*)tkrDigiIter.Next()) ) {
		
		int Tower = pTkrDigi->getTower().id();
		int Layer = pTkrDigi->getBilayer();
		int View  = pTkrDigi->getView();

		int Plane = 2*Layer+View;
		
		if ( DEBUG ) std::cout << "tower, bilayer, view, plane: " <<Tower<<" "<<Layer<<" "<<View<<" "<<Plane<<std::endl;
		const int TOT[2] = { pTkrDigi->getToT(0), pTkrDigi->getToT(1) };
		
		TrueToT0[Plane] = TOT[0];
		TrueToT1[Plane] = TOT[1];

		if ( DEBUG ) std::cout << "ToT0 = : " << TrueToT0[Plane] << " ToT1 = " << TrueToT1[Plane] << std::endl;
		//last0strip = pTkrDigi->getLastController0Strip();
		//		if ( DEBUG ) std::cout << "last strip of controller 0: " << last0strip << std::endl;
		
		TrueTkrNumHits[Plane] = pTkrDigi->getNumHits();
		TkrTotalNumHits += TrueTkrNumHits[Plane];
		
		if ( DEBUG ) std::cout << "num hits: " << TrueTkrNumHits[Plane]  << std::endl;
		
		if ( TrueTkrNumHits[Plane] > 0 ) 
		  {
		    for ( int h=0; h<TrueTkrNumHits[Plane] ; h++ ) 
		      TrueTkrHits[Plane][h] = pTkrDigi->getHit(h);
		    		    
		    if ( DEBUG ) 
		      { 
			std::cout << "list of hits:";
			for ( int h=0; h<TrueTkrNumHits[Plane] ; h++ ) std::cout << " " << TrueTkrHits[Plane][h];
			std::cout << std::endl;
		      }
		    
		  }
		
	      }
	    }
	  
	  
	  //////////////////////////////////////////////////
	}
      // RECON ONLY analysis
      if (rec) {  // if we have recon data proccess it
	TkrRecon *tkrRec = rec->getTkrRecon();
	if (tkrRec)
	  {
	    ////////////////// TKR CLUSTER: ////////////////////////////////
	    const TObjArray *clusCol = tkrRec->getClusterCol();
	    TkrNumClus = clusCol->GetEntries();
	    
	    TIter tkrClusIter(clusCol);
	    TkrCluster* pTkrClus = 0;
	    int clusIdx=0;
	    while ( (pTkrClus=(TkrCluster*)tkrClusIter.Next()) ) 
	      {
		TkrClusX[clusIdx]     = pTkrClus->getPosition().X();
		TkrClusY[clusIdx]     = pTkrClus->getPosition().Y();
		TkrClusZ[clusIdx]     = pTkrClus->getPosition().Z();
		TkrClusPlane[clusIdx] = pTkrClus->getPlane();
		TkrClusView[clusIdx]  = pTkrClus->getView();
		clusIdx++;
	      }
	    ////////////////// TKR TRACKS: ////////////////////////////////
	    const TObjArray *tkrCol = tkrRec->getTrackCol();
	    TkrNumTracks = tkrCol->GetEntries();
	    
	    if ( TkrNumTracks > 0 ) {
	      const TkrKalFitTrack* track1 = (TkrKalFitTrack*)tkrCol->First();
	      const int nHits = track1->getNumHits();
	      for ( int i=0; i<nHits; ++i )
		TkrTrk1Clusters[i] = track1->getHitPlane(i)->getIdHit();
	    }
	    
	    ////////////////// TLKR VERTEX: ////////////////////////////////
	    
	    const TObjArray *vtxCol = tkrRec->getVertexCol();
	    
	    TkrNumVtx = vtxCol->GetEntries();
	    
	    TkrVertex *vtx = (TkrVertex*) vtxCol->First();
	    
	    if(vtx)
	      {
		TVector3 dir = vtx->getDirection();
		
		TkrVtxX = dir.X();
		TkrVtxY = dir.Y();
		TkrVtxZ = dir.Z();
		Theta = dir.Theta();
		Phi =dir.Phi();
		if(dir.Z()==0)
		  {
		    
		    ThetaXZ = 0.0;
		    ThetaYZ = 0.0;
		  }
		else
		  {
		    ThetaXZ = dir.X()/dir.Z();
		    ThetaYZ = dir.Y()/dir.Z();
		  }
		
		if ( DEBUG ) std::cout<<TkrVtxX<<" "<<TkrVtxY<<" "<<TkrVtxZ<<std::endl;
	      }// vertex
	  } // tkrrecon
      } //recon data
      if ( DEBUG ) std::cout << "num tkrDigi before filling..." << nTkrDigi << std::endl;

      TIter TreeIter(TreeCollection);
      TTree* aTree = 0;
      int i=0;
      
      while ( (aTree=(TTree*)TreeIter.Next()) )
	{
	  TkrNumHits = TrueTkrNumHits[i];
	  ToT0       = TrueToT0[i];
	  ToT1       = TrueToT1[i];
	  //	  TkrHits    = TrueTkrHits[i];
	  if(DIAGN)
	    {
	      const char* name = aTree->GetName();
	      const int index0 = getIndex(name, false);
	      const int index1 = getIndex(name, true);
	      TriggerReq0 = TkrDiagnostics[index0];
	      TriggerReq1 = TkrDiagnostics[index1];
	      if ( DEBUG ) std::cout << "loop " << std::setw(8) << name << std::setw(3) << index0 << std::setw(3) << index1 << ' '
				     << TriggerReq0 << ' ' << TriggerReq1 << ' ' << TkrNumHits << std::endl;
	    }
	  if (TkrNumHits > 0 || i==0 ) 
	    {
	      for ( int h=0; h<TkrNumHits; h++ ) 
		TkrHits[h] = TrueTkrHits[i][h];
	      if ( DEBUG ) std::cout<<"fill "<<i<<" "<<TkrNumHits<<std::endl;
	    }
	  aTree->Fill();
	  i++;
	  
	}
      //      if(nTkrDigi>0) myTree->Fill();
      if ( DEBUG ) std::cout << "...filled!" << nTkrDigi << std::endl;
    }  // end analysis code in event loop
  
  m_StartEvent = curI;
  if (m_TreeFileName=="") m_TreeFileName="MyRootFile.root";
  
  TFile myFile("MyRootFile.root","RECREATE","TreeFile",1);
  TreeCollection->Write();
  myFile.Close();
  std::cout<<m_TreeFileName<<std::endl;
}

#endif
  
