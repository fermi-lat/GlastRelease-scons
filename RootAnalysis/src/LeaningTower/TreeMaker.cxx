#include "TreeMaker.h"

#include <vector>
#include <bitset>
#include <iomanip>
#include <map>
#include <utility>

const int MaxNumPlanes = 36;
const Int_t NumGTCC = 8;
const Int_t NumGTRC = 9;
const Int_t MaxNumTkrHits = 128;
const Int_t MaxNumTkrClusters = 128;

/*
 * The order of read-out of the GTCC's (LAT-TD-00605 p.106
 * http://www-glast.slac.stanford.edu/IntegrationTest/ONLINE/docs/TEM.pdf):
 * 6, 3, 7, 2, 5, 0, 4, 1
 */
const int positionOfGTCC[8] /*in the TCloneArray*/ = { 5, 7, 3, 1, 6, 4, 0, 2 };

const bool DEBUG = 0;

#define DIAGN

// returns position of GTCC i in the TCloneArray m_TkrDiagnostics of the digi
int indexToGTCC ( int i ) { return positionOfGTCC[i]; }

// returns the index of a GTCC/GTRC pair in the TkrDiagnostics array
int getIndex ( int GTCC, int GTRC ) { return GTCC*NumGTRC+GTRC; }

// returns index of a plane/ROC pair in TkrDiagnostics. Left is 0, right is 1.
int getIndex ( const TString planeName, const bool side ) {
    // assume e.g. "LayerX12"
    static const int first = 5;
    char view = planeName[first];
    int num = atoi(planeName(first+1,planeName.Length()-first-1).Data());
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
    
    return getIndex(GTCC, GTRC);
}

void TreeMaker::CreateTree(Int_t numEvents) {
    // variables to be stored in "Header"
    Int_t EventId;
    Int_t RunId;
    Int_t TemId;
    Int_t TkrTotalNumHits;
    Double_t EbfTime;  // here double might be needed
    Int_t LevelOneTrigger;
    Int_t TkrDiagnostics[NumGTCC*NumGTRC]; // this could be sorted into planes

    // variables to be stored in every "Layer"
    Int_t ToT0;
    Int_t ToT1;
    Bool_t TriggerReq0;
    Bool_t TriggerReq1;
    Int_t TkrNumHits;
    Int_t TkrHits[MaxNumTkrHits];  // strip id's
    // temporary store for later filling of the Layers
    Int_t ToT0Array[MaxNumPlanes];
    Int_t ToT1Array[MaxNumPlanes];
    Int_t TkrNumHitsArray[MaxNumPlanes];
    Int_t TkrHitsArray[MaxNumPlanes][MaxNumTkrHits];

    // variables to be stored in Recon
    // clusters
    Int_t TkrNumClus;
    Int_t TkrClusLayer[MaxNumTkrClusters]; // this is really layer, not plane
    Int_t TkrClusView[MaxNumTkrClusters];
    Float_t TkrClusX[MaxNumTkrClusters];
    Float_t TkrClusY[MaxNumTkrClusters];
    Float_t TkrClusZ[MaxNumTkrClusters];
    //  tracks
    Int_t TkrNumTracks;
    Int_t TkrTrk1NumClus;
    Int_t TkrTrk1Clusters[MaxNumTkrClusters];
    // vertex
    Int_t TkrNumVtx;
    Float_t TkrVtx1X, TkrVtx1Y, TkrVtx1Z;
    Float_t TkrVtx1Theta, TkrVtx1Phi, TkrVtx1ThetaXZ, TkrVtx1ThetaYZ;

    UInt_t TkrDigi3RowBits,TkrTrgReq3RowBits;

    TFile myFile(m_TreeFileName, "RECREATE", "TreeFile", 1);
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
	digiTree->SetBranchStatus("m_gem", 1);
#ifdef DIAGN
        digiTree->SetBranchStatus("m_levelOneTrigger", 1);
        digiTree->SetBranchStatus("m_ebfTime*", 1);
        digiTree->SetBranchStatus("m_numTkrDiagnostics", 1);
#endif
        for(int i=0;i<MaxNumPlanes;i++) {
            TString TS = "Layer";
            TS += i%2 == 0 ? "X" : "Y";
            TS += i/2;
            if ( DEBUG ) std::cout << TS << std::endl;
            TTree *Layer = new TTree(TS,TS);
            Layer->Branch("ToT0",&ToT0,"ToT0/I");
            Layer->Branch("ToT1",&ToT1,"ToT1/I");
#ifdef DIAGN
            Layer->Branch("TriggerReq0",&TriggerReq0,"TriggerReq0/B");
            Layer->Branch("TriggerReq1",&TriggerReq1,"TriggerReq1/B");
#endif
            Layer->Branch("TkrNumHits",&TkrNumHits,"TkrNumHits/I");
            Layer->Branch("TkrHits",TkrHits,"TkrHits[TkrNumHits]/I");
            TreeCollection->Add((TObject*)Layer);
        }
        TTree *Header = new TTree("Header","Header");
        Header->Branch("EventId",&EventId,"EventId/I");
        Header->Branch("RunId",&RunId,"RunId/I");
        Header->Branch("TkrTotalNumHits",&TkrTotalNumHits,"TkrTotalNumHits/I");
        Header->Branch("TemId",&TemId,"TemId/I");
#ifdef DIAGN
        Header->Branch("EbfTime",&EbfTime,"EbfTime/D");
        Header->Branch("LevelOneTrigger",&LevelOneTrigger,"LevelOneTrigger/I");
        TString tag("TkrDiagnostics[");
        tag += NumGTCC * NumGTRC;
        tag += "]/I";
        Header->Branch("TkrDiagnostics", TkrDiagnostics, tag);
	Header->Branch("TkrTrgReq3RowBits",&TkrTrgReq3RowBits,"TkrTrgReq3RowBits/i");
	Header->Branch("TkrDigi3RowBits",&TkrDigi3RowBits,"TkrDigi3RowBits/i");
#endif
        TreeCollection->Add(Header);
    }
  
    if ( reconTree ) {
        reconTree->SetBranchStatus("*",0);  // disable all branches
        // activate desired branches
        reconTree->SetBranchStatus("m_cal*", 1);  
        reconTree->SetBranchStatus("m_tkr*", 1);
        reconTree->SetBranchStatus("m_acd*", 1);
        reconTree->SetBranchStatus("m_eventId", 1); 
        reconTree->SetBranchStatus("m_runId", 1);
    
        TTree* Recon = new TTree("Recon", "Recon");
        ////////// Clusters: ////////////////////////////////////////
        Recon->Branch("TkrNumClus",&TkrNumClus,"TkrNumClus/I");
        Recon->Branch("TkrClusX",TkrClusX,"TkrClusX[TkrNumClus]/F");
        Recon->Branch("TkrClusY",TkrClusY,"TkrClusY[TkrNumClus]/F");
        Recon->Branch("TkrClusZ",TkrClusZ,"TkrClusZ[TkrNumClus]/F");
        Recon->Branch("TkrClusLayer",TkrClusLayer,"TkrClusLayer[TkrNumClus]/I");
        Recon->Branch("TkrClusView",TkrClusView,"TkrClusView[TkrNumClus]/I");
        /////////// Tracks: ////////////////////////////////////////
        Recon->Branch("TkrNumTracks",&TkrNumTracks,"TkrNumTracks/I");
        Recon->Branch("TkrTrk1NumClus",&TkrTrk1NumClus,"TkrTrk1NumClus/I");
        Recon->Branch("TkrTrk1Clusters",TkrTrk1Clusters,
                      "TkrTrk1Clusters[TkrNumClus]/I");
        /////////// Vertex: ////////////////////////////////////////
        Recon->Branch("TkrNumVtx",&TkrNumVtx,"TkrNumVtx/I");
        Recon->Branch("TkrVtx1X",&TkrVtx1X,"TkrVtx1X/F");
        Recon->Branch("TkrVtx1Y",&TkrVtx1Y,"TkrVtx1Y/F");
        Recon->Branch("TkrVtx1Z",&TkrVtx1Z,"TkrVtx1Z/F");
        Recon->Branch("TkrVtx1Theta",&TkrVtx1Theta,"TkrVtx1Theta/F");
        Recon->Branch("TkrVtx1Phi",&TkrVtx1Phi,"TkrVtx1Phi/F");
        Recon->Branch("TkrVtx1ThetaXZ",&TkrVtx1ThetaXZ,"TkrVtx1ThetaXZ/F");
        Recon->Branch("TkrVtx1ThetaYZ",&TkrVtx1ThetaYZ,"TkrVtx1ThetaYZ/F");
        TreeCollection->Add(Recon);
    }

    // determine how many events to process
    Int_t numEntries = GetEntries();
    std::cout << "\nNum Events in File is: " << numEntries << std::endl;
    numEvents = std::min(numEvents, numEntries);

    //////////////////////////////////////////////////
    // BEGINNING OF EVENT LOOP

    // I encountered a case where the digi file made from ldf contains events
    // which don't have a 3-in-a-row trigger.  In that case, the recon file
    // made from Gleam misses events.  Thus, I have to introduce an event record
    // counter per tree.
    UInt_t evtCounterMc    = 0;
    UInt_t evtCounterDigi  = 0;
    UInt_t evtCounterRecon = 0;
    Int_t digiEventIdOld, reconEventIdOld, mcEventIdOld;
    digiEventIdOld = reconEventIdOld = mcEventIdOld = -1;
    for ( Int_t evtCounter=0; evtCounter<numEvents; ++evtCounter ) {
        if(DEBUG) std::cout << "event loop counter: " << evtCounter <<std::endl;
        // Reset some arrays and variables
        int nTkrDigi = 0;
        TkrTotalNumHits = 0;
        for ( int i=0; i<MaxNumPlanes; i++) {
            ToT0Array[i] = ToT1Array[i] = -1;
            TkrNumHitsArray[i] = 0;
        }
        TkrNumClus = 0;
        TkrNumTracks = 0;
        TkrTrk1NumClus = 0;

        // some progress bar
        int fraction = (int)(100.0 * evtCounter / numEvents);
        static int target_fraction = 0;
        if ( fraction > target_fraction ) {
            target_fraction = fraction;
            if ( fraction >= 10 )
                target_fraction = target_fraction/10*10 + 9;
            std::cout << fraction <<"% complete: event "<<evtCounter<<std::endl;
        }
      
        if (mc)
            mc->Clear();
        if (evt)
            evt->Clear();
        if (rec)
            rec->Clear();
      
        Int_t digiEventId, reconEventId, mcEventId;
        digiEventId = reconEventId = mcEventId = 100000000;
        UInt_t digiRunNum, reconRunNum, mcRunNum;
        digiRunNum = reconRunNum = mcRunNum = 0;
        GetEvent(evtCounterMc, evtCounterDigi, evtCounterRecon);


	//Giving up on this event if more than one tower is involved.
	Int_t iTem;
	Int_t iCount=0;
	Gem gem = evt->getGem();
        UShort_t tkrVec = gem.getTkrVector();
        for (iTem = 0; iTem < 16; iTem++) {
	  Int_t offset = ((tkrVec>>iTem)&1);
	  //	  std::cout << "Tem "<< iTem << " " << offset  << std::endl;
	  if (offset==1){
	    iCount++;
	    TemId = iTem;
	  }
        }
	//flagging events for which more than one TEM contribute : to be discarded when TreeMaker ROOT file is read.
	//For some reason putting a continue here screws up the rest of the event loop (same event repeated until the end)
       	if(iCount>1) 
	  TemId=-1;

	  
        // synchronizing the three trees
        if (mc)
            mcEventId = mc->getEventId();
        if (evt)
            digiEventId = evt->getEventId();
            if (rec)
            reconEventId = rec->getEventId();
        // agreeing on a common event id
        EventId = std::min(std::min(mcEventId, digiEventId), reconEventId);

        /*
        std::cout << "event ids: " << mcEventId << ' ' << digiEventId << ' '
                  << reconEventId << ' ' << EventId << std::endl;
        std::cout << "old event ids: " << mcEventIdOld << ' ' << digiEventIdOld
                  << ' ' << reconEventIdOld << std::endl;
        */

        //////////////////////////////////////////////////
        // Monte Carlo tree
        if ( mc && mcEventId != EventId )
            std::cout << "skipping mc tree: " << EventId << ' ' << mcEventId
                      << ' ' << digiEventId << ' ' << reconEventId << std::endl;
        else if ( mc ) {  // if we have mc data process it
            ++evtCounterMc;
            mcRunNum  = mc->getRunId();
            if ( DEBUG ) std::cout << "MC: run number: " << mcRunNum
                                   << "  event id: " << mcEventId << std::endl;
        }

        // Digi tree
        if ( evt && digiEventId != EventId )
            std::cout << "skipping digi tree: " << EventId << ' ' << mcEventId
                      << ' ' << digiEventId << ' ' << reconEventId << std::endl;
        else if ( evt ) {  // if we have digi data process it

            ++evtCounterDigi;
            digiRunNum = evt->getRunId();
            if ( DEBUG ) std::cout << "DIGI: run number: " << digiRunNum
                                   << "  event id: " << digiEventId <<std::endl;
            RunId = digiRunNum;
#ifdef DIAGN
            // we shouldn't do this if the digi is generated from mc
            L1T l1t = evt->getL1T();
	    if ( DEBUG )
                l1t.Print();
            LevelOneTrigger = l1t.getTriggerWord() & 0x1f;
            // there are more bits set, but why?
            bool tkrTrigger = l1t.getTkr3InARow();
	    TkrTrgReq3RowBits = evt->getL1T().getTrgReqTriRowBits(0);
	    TkrDigi3RowBits = evt->getL1T().getDigiTriRowBits(0);
            if ( DEBUG )
                std::cout << "bool 3-in-a-row: " << tkrTrigger
                          << "  TkrTrgReq3RowBits: "
                          << std::setw(6) << TkrTrgReq3RowBits << ' '
                          << static_cast<std::bitset<sizeof(UInt_t)*8> >
                    (TkrTrgReq3RowBits)
                          << "  TkrDigi3RowBits: "
                          << std::setw(6) << TkrDigi3RowBits << ' '
                          << static_cast<std::bitset<sizeof(UInt_t)*8> >
                    (TkrDigi3RowBits)
                          <<std::endl;

            // Ebf time
            static Double_t EbfTimeStart = evt->getEbfTimeSec();
            EbfTime = ( evt->getEbfTimeSec() - EbfTimeStart )
                + evt->getEbfTimeNanoSec() * 1E-9;
            if ( DEBUG ) std::cout << "EbfTime: " << EbfTime << std::endl;

            // and now, for the TkrDiagnostics
            int numGTCCentries = evt->getTkrDiagnosticCol()->GetEntries();
            if ( NumGTCC != numGTCCentries ) {
                static bool NumGTCCwarning = true;
                if ( DEBUG || NumGTCCwarning ) // print at least once
                    // e.g., for Monte Carlo data
                    std::cerr << "NumGTCC=" << NumGTCC
                              << " != evt->getTkrDiagnosticCol()->GetEntries()="
                              << numGTCCentries << std::endl;
                NumGTCCwarning = false;
            }
            else
                for ( int GTCC=0; GTCC<NumGTCC; ++GTCC ) { 
                    std::bitset<NumGTRC> word = 0;
                    if ( tkrTrigger && GTCC < numGTCCentries )
                        word =
                        evt->getTkrDiagnostic(indexToGTCC(GTCC))->getDataWord();
                    if ( DEBUG ) std::cout << "TkrDiagnosticData[" << GTCC
                                           << "] word(" << word << ") ";
                    for ( int GTRC=NumGTRC-1; GTRC>=0; --GTRC ) {
                        TkrDiagnostics[getIndex(GTCC,GTRC)] = word[GTRC];
                        if ( DEBUG )
                         std::cout<<(TkrDiagnostics[getIndex(GTCC,GTRC)]>0)?1:0;
                    }
                    if ( DEBUG ) std::cout << std::endl;
                }
#endif
            // TkrDigi
            const TObjArray* tkrDigiCol = 0;
            tkrDigiCol = evt->getTkrDigiCol();
            //////////////////////////////////////////////////
            if ( tkrDigiCol ) {
                nTkrDigi = tkrDigiCol->GetEntries();
                if (DEBUG) std::cout << "num tkrDigi " << nTkrDigi << std::endl;

                TIter tkrDigiIter(tkrDigiCol);
                TkrDigi* pTkrDigi = 0;
                while ( (pTkrDigi=(TkrDigi*)tkrDigiIter.Next()) ) {
                    int Tower = pTkrDigi->getTower().id();
                    int Layer = pTkrDigi->getBilayer();
                    int View  = pTkrDigi->getView();
                    // this is not a "Ritz" plane, so I call it Nicolas plane
                    int NicolasPlane = 2*Layer+View;
                    if ( DEBUG ) std::cout << "tower, bilayer, view, plane: "
                                           <<Tower<<" "<<Layer<<" "<<View<<" "
                                           <<NicolasPlane<<std::endl;
                    ToT0Array[NicolasPlane] = pTkrDigi->getToT(0);
                    ToT1Array[NicolasPlane] = pTkrDigi->getToT(1);
                    if ( DEBUG ) std::cout << "ToT0 = "<<ToT0Array[NicolasPlane]
                          << " ToT1 = " << ToT1Array[NicolasPlane] << std::endl;
                    //last0strip = pTkrDigi->getLastController0Strip();

                    TkrNumHitsArray[NicolasPlane] = pTkrDigi->getNumHits();
                    if ( DEBUG ) std::cout << "num hits: "
                                  << TkrNumHitsArray[NicolasPlane] << std::endl;
                    TkrTotalNumHits += TkrNumHitsArray[NicolasPlane];
		
                    if ( TkrNumHitsArray[NicolasPlane] > 0 ) {
                        // filling the array of TkrHits.
                        // Don't fill more than MaxNumTkrHits!
                        for ( int h=0; h<TkrNumHitsArray[NicolasPlane]; h++ ) {
                            if ( h >= MaxNumTkrHits ) {
                                std::cerr << "Event " << EventId << " has "
                                          << TkrNumHitsArray[NicolasPlane]
                                          << " hits in layer " << Layer
                                          << " view " << View << std::endl;
                                break;
                            }
                            TkrHitsArray[NicolasPlane][h] = pTkrDigi->getHit(h);
                        }
                        if ( DEBUG ) { 
                            std::cout << "list of hits:";
                            for ( int h=0; h<TkrNumHitsArray[NicolasPlane]; h++)
                                std::cout<<' '<<TkrHitsArray[NicolasPlane][h];
                            std::cout << std::endl;
                        }
                    }
                }
            }
        }

        //////////////////////////////////////////////////
        // RECON tree
        if ( rec && reconEventId != EventId )
            std::cout << "skipping recon tree: " << EventId << ' ' << mcEventId
                      << ' ' << digiEventId << ' ' << reconEventId << std::endl;
        else if ( rec ) {  // if we have recon data proccess it
	    ++evtCounterRecon;
            reconRunNum = rec->getRunId();
            TkrRecon* tkrRec = rec->getTkrRecon();
            if ( tkrRec ) {
                ////////////////// TKR TRACKS: ////////////////////////////////
                const TObjArray *tkrCol = tkrRec->getTrackCol();
                TkrNumTracks = tkrCol->GetEntries();
		const TkrTrack* track1 = 0;
		int flag = 1;
		int count = 0;
                if ( TkrNumTracks > 0 ) {
		  track1 = (TkrTrack*)tkrCol->First();
		  TkrTrk1NumClus = track1->GetEntries();
 		}
		////////////////// TKR CLUSTER: ////////////////////////////////
                const TObjArray* clusCol = tkrRec->getClusterCol();
                TkrNumClus = clusCol->GetEntries();
                TIter tkrClusIter(clusCol);
                TkrCluster* pTkrClus = 0;
                int clusIdx = 0;
		int trk1Idx = 0;
                // what was eps for?  Comment now, but to be removed one day!
                //		float eps = 0.0001;
                while ( ( pTkrClus = (TkrCluster*)tkrClusIter.Next() ) ) {
                    // filling the arrays of TkrClus.
                    // Don't fill more than MaxNumTkrClusters!
                    if ( clusIdx >= MaxNumTkrClusters ) {
                        std::cerr << "Event " << EventId << " has "
                                  << TkrNumClus << " clusters" << std::endl;
                        TkrNumClus = MaxNumTkrClusters;
                        break;
                    }
                    TkrClusX[clusIdx]     = pTkrClus->getPosition().X();
                    TkrClusY[clusIdx]     = pTkrClus->getPosition().Y();
                    TkrClusZ[clusIdx]     = pTkrClus->getPosition().Z();
                    // getPlane()! doesn't return the plane but the
                    // layer index (number 0 - 17)
                    TkrClusLayer[clusIdx] = pTkrClus->getLayer();
                    TkrClusView[clusIdx]  = pTkrClus->getTkrId().getView();

		    //now map the cluster to the best track:
		    if( TkrNumTracks > 0){
		      TIter trk1HitsItr(track1);
		      TkrTrackHit* pTrk1Hit = 0;
		      int nCount =-1;
		      while( (pTrk1Hit = (TkrTrackHit*)trk1HitsItr.Next()) )
			{
			  nCount++;
			  const TkrCluster* pTrk1Clus = pTrk1Hit->getClusterPtr();
			  if(pTrk1Clus)
			    {
			      //                              std::cout<<pTrk1Clus->getPosition().X() - pTkrClus->getPosition().X()<<std::endl;
			      //                              std::cout<<pTrk1Clus->getPosition().Y() - pTkrClus->getPosition().Y()<<std::endl;
			      //                              std::cout<<pTrk1Clus->getPosition().Z() - pTkrClus->getPosition().Z()<<std::endl;
			      //                              std::cout<<pTrk1Clus->getLayer()        - TkrClusLayer[clusIdx]<<std::endl;
			      //			      std::cout<<pTrk1Clus->getTkrId().getView() -  TkrClusView[clusIdx]<<std::endl;
			      if(
				 ( pTrk1Clus->getPosition().X()       ==  pTkrClus->getPosition().X() )
				 && ( pTrk1Clus->getPosition().Y()    ==  pTkrClus->getPosition().Y() )
				 && ( pTrk1Clus->getPosition().Z()    ==  pTkrClus->getPosition().Z() )
				 && ( pTrk1Clus->getLayer()           == pTkrClus->getLayer() )
				 && ( pTrk1Clus->getTkrId().getView() == pTkrClus->getTkrId().getView() )
				 )
				{
				  TkrTrk1Clusters[trk1Idx] = clusIdx;
				  //                                  std::cout<<"track hit "<< nCount <<" matches cluster " <<clusIdx <<std::endl;
				  trk1Idx++;
				  break;
				}
			      else
 				{
				  continue;
				  // 				  std::cout<<"cluster "<<clusIdx <<" does not match track hit " <<nCount<<std::endl;
 				}
			    }
 			  else
 			    {
			      //			      std::cout<<"Track hit "<<nCount << "is a gap"<<std::endl;
			      if(flag==1)count++;
			      // 			      TkrTrk1Clusters[trk1Idx] = -1;
			      // 			      trk1Idx++;
 			    }
			}
		    }
		    flag = 0;
                    ++clusIdx;    
                }
		TkrTrk1NumClus -=count;
                ////////////////// TKR TRACKS: ////////////////////////////////
//                 const TObjArray *tkrCol = tkrRec->getTrackCol();
//                 TkrNumTracks = tkrCol->GetEntries();

//                 if ( TkrNumTracks > 0 ) {
//                     const TkrTrack* track1 = (TkrTrack*)tkrCol->First();
//                     TkrTrk1NumClus = track1->Size();
//                     for ( int i=0; i<TkrTrk1NumClus; ++i ) {
//                         if ( i >= MaxNumTkrClusters ) {
//                             std::cerr << "Event " << EventId
//                                       << " has a first track with "
//                                       << TkrTrk1NumClus<<" clusters"<<std::endl;
//                             break;
//                         }
// 			TIter trk1HitsItr(track1);
// 			TkrTrackHit* pTrk1Hit = 0;
// 			while( (pTrk1Hit = (TkrTrackHit*)trk1HitsIter.Next()) )
// 			  {
			    
// 			  }
// 			//                        TkrTrk1Clusters[i] = track1->getHitPlane(i)->getIdHit();
//                     }
//                 }
		

                ////////////////// TKR VERTEX: ////////////////////////////////
                const TObjArray *vtxCol = tkrRec->getVertexCol();
                TkrNumVtx = vtxCol->GetEntries();
                TkrVertex *vtx = (TkrVertex*) vtxCol->First();
	    
                if ( vtx ) {
                    TVector3 dir = vtx->getDirection();
                    TkrVtx1X = dir.X();
                    TkrVtx1Y = dir.Y();
                    TkrVtx1Z = dir.Z();
                    TkrVtx1Theta = dir.Theta();
                    TkrVtx1Phi =dir.Phi();
                    if(dir.Z()==0) {
                        TkrVtx1ThetaXZ = 0.0;
                        TkrVtx1ThetaYZ = 0.0;
                    }
                    else {
                        TkrVtx1ThetaXZ = dir.X()/dir.Z();
                        TkrVtx1ThetaYZ = dir.Y()/dir.Z();
                    }
                    if ( DEBUG ) std::cout<<TkrVtx1X<<" "<<TkrVtx1Y<<" "
                                          <<TkrVtx1Z<<std::endl;
                }// vertex
            } // tkrrecon
        } //recon data

        if ( DEBUG ) std::cout << "num tkrDigi before filling..." << nTkrDigi
                               << std::endl;

        // now we are filling all trees

        /*
        static bool flag = false;
        if ( digiEventId > 131060 )
            flag = true;
        if ( flag ) {
            std::cout << "event ids: " << mcEventId << ' ' << digiEventId << ' ' << reconEventId << ' ' << EventId << std::endl;
            std::cout << "old event ids: " << mcEventIdOld << ' ' << digiEventIdOld << ' ' << reconEventIdOld << std::endl;
            std::cout << "run ids: " << mcRunNum << ' ' << digiRunNum << ' ' << reconRunNum << ' ' << RunId << std::endl;
        }
        */
        TIter TreeIter(TreeCollection);
        TTree* aTree = 0;
        int iLayer = 0;
        while ( ( aTree=(TTree*)TreeIter.Next() ) ) {
            // it would be saver to use std::map.
            // If nothing changes, the sequence of planes is:
            // 0-35: LayerX0 - LayerY17
            // 36:   Header
            // 37:   Recon
            if ( iLayer < 36 ) { // the tree is a Layer
                TkrNumHits = TkrNumHitsArray[iLayer];
                ToT0       = ToT0Array[iLayer];
                ToT1       = ToT1Array[iLayer];
                //	  TkrHits    = TkrHitsArray[iLayer];
                // if (TkrNumHits>0 || iLayer==0 ) {//what serves the iLayer==0?
                if ( TkrNumHits > 0 ) {
                    // copying up to MaxNumTkrHits TkrHits
                    for ( int h=0; h<std::min(MaxNumTkrHits,TkrNumHits); h++ ) 
                        TkrHits[h] = TkrHitsArray[iLayer][h];
                    if ( DEBUG ) std::cout << "fill " << iLayer << ' '
                                           << TkrNumHits << std::endl;
                }
#ifdef DIAGN
                const char* name = aTree->GetName();
                const int index0 = getIndex(name, false);
                const int index1 = getIndex(name, true);
                TriggerReq0 = TkrDiagnostics[index0];
                TriggerReq1 = TkrDiagnostics[index1];
                if ( DEBUG ) std::cout << "loop " << std::setw(8) << name
                                       << std::setw(3) << index0
                                       << std::setw(3) << index1 << ' '
                                       << TriggerReq0 << ' ' << TriggerReq1
                                       << ' ' << TkrNumHits << std::endl;
#endif
            }
            aTree->Fill();
            ++iLayer;
        }

        if ( DEBUG ) std::cout << "...filled!" << nTkrDigi << std::endl;
	
    }  // end analysis code in event loop

    TreeCollection->Write(0, TObject::kWriteDelete);
    myFile.Close();
    std::cout << "m_TreeFileName " << m_TreeFileName << std::endl;
}
