// -*- Mode: c++ -*-
#ifndef CelUtil_cxx
#define CelUtil_cxx
/*
* Project: GLAST
* Package: rootUtil
*    File: $Id$
* Authors:
*   EC, Eric Charles,    SLAC              echarles@slac.stanford.edu
*
* Copyright (c) 2007
*                   Regents of Stanford University. All rights reserved.
*
*
*
*/

// This Class's header
#include "rootUtil/CelUtil.h"

// c++/stl headers
#include <iostream>

// ROOT Headers
#include <TTree.h>
#include <TChain.h>
#include <TFile.h>
#include <TObjArray.h>

// Other headers from this package
#include "rootUtil/CelEventComponent.h"
#include "rootUtil/CompositeEventList.h"
#include "rootUtil/FileUtil.h"


CompositeEventList* CelUtil::mergeCelFiles(TCollection& skimFiles, const char* fileName, const char* option) {
  // Merge a Collection of input files into a single composite event list

  TList cels;
  TIterator* itr = skimFiles.MakeIterator();
  TObject* aFileObj(0);
  while ( (aFileObj = itr->Next()) != 0 ) {
    TFile* aFile = dynamic_cast<TFile*>(aFileObj);
    if ( aFile == 0 ) {
      std::cerr << "Input object not a file " << aFileObj->GetName() << std::endl;
      return 0;
    }

    // Open the skim
    CompositeEventList* nextSkim = new CompositeEventList;
    TFile* check = nextSkim->openCelFile(aFile->GetName());
    if ( check == 0 ) {
      std::cerr << "Failed to open skim file " << aFile->GetName() << std::endl;
      return 0;
    }
    // Append it to the list of cels
    cels.AddLast(nextSkim);
  }
  // Call the other mergre method
  return mergeCompositeEventLists(cels,fileName,option);
}


// Merge a Collection of composite event lists
CompositeEventList* CelUtil::mergeCompositeEventLists
 ( TCollection & cels, const char* fileName, const char* option)
 {

  TList entryTreeList ;
  TList fileTreeList ;

  // Create the output File if requested
  TFile* outFile(0);
  if ( fileName != 0 ) {
    outFile = TFile::Open(fileName,option);
    if ( 0 == outFile ) {
      std::cerr << "Failed to open output file " << fileName << std::endl;
      return 0;
    } else {
      std::cout << "Opened output file " << outFile->GetName() << std::endl;
    }
  }

  // The Event Tree has to be re-done
  /// REWRITE THIS WITH CelEventLink
  Long64_t* eventIndex = new Long64_t(-1);
  Long64_t* fileSetIndexOut  = new Long64_t(-1);
  Long64_t* fileSetOffset = new Long64_t(0);
  TTree* linkTreeOut = new TTree("EventLinks","EventLinks");
  linkTreeOut->Branch("Event_Index",(void*)eventIndex,"Event_Index/L");
  linkTreeOut->Branch("Event_SetIndex",(void*)fileSetIndexOut,"Event_SetIndex/L");
  linkTreeOut->Branch("Event_SetOffset",(void*)fileSetOffset,"Event_SetOffset/L");

  // Loop on cels
  TIterator* itr = cels.MakeIterator();
  TObject* aCelObj(0);
  while ( (aCelObj = itr->Next()) != 0 ) {
    CompositeEventList* aCel = dynamic_cast<CompositeEventList*>(aCelObj);
    if ( 0 == aCel ) {
      std::cerr << "Input object not a CompositeEventList " << aCelObj->GetName() << std::endl;
      return 0;
    }    
    
    // Add the event & file trees to their respective lists
    entryTreeList.AddLast(aCel->entryTree());
    fileTreeList.AddLast(aCel->fileTree());    

    Long64_t fileSetIndexIn(0);
    Long64_t fileSetIndexSave(-1);
    TTree* linkTree = aCel->linkTree();
    linkTree->SetBranchAddress("Event_SetIndex",&fileSetIndexIn);

    // Redo the data in the Event tree
    Long64_t nEvtCurrent = linkTree->GetEntries();
    for ( Long64_t iEvt(0); iEvt < nEvtCurrent; iEvt++ ) {
      Int_t nB = linkTree->GetEntry(iEvt);
      if ( nB < 0 ) {
	// There was a problem, clean up
	std::cerr << "Failed to read entry in input link tree " << std::endl;
	if ( outFile != 0 ) {
	  outFile->Close();
	  delete outFile;
	}
	delete eventIndex;
	delete fileSetIndexOut;
	delete fileSetOffset;
	delete itr;
	delete linkTree;
	return 0;
      }
      // Event Counter goes up by one
      *eventIndex += 1;
      if ( fileSetIndexIn != fileSetIndexSave ) {
	// new meta data entry, latch values
	fileSetIndexSave = fileSetIndexIn;
	*fileSetIndexOut += 1;
	*fileSetOffset = *eventIndex;
      }
      // Fill this entry in the link tree
      linkTreeOut->Fill();
    }

  }
  delete itr;
  delete eventIndex;
  delete fileSetIndexOut;
  delete fileSetOffset;

  // Just Merge the Event and File trees the Normal way
  TTree* entryTree = TTree::MergeTrees(&entryTreeList);
  TTree* fileTree = TTree::MergeTrees(&fileTreeList);

  // Build the Pointer Skim Object
  CompositeEventList * theCel = new CompositeEventList(*entryTree,*linkTreeOut,*fileTree) ;

  // Write the Merged File, if requested
  if ( outFile != 0 )
   {
    outFile->Write();
    outFile->Close();
   }
  return theCel ;  
 }


#endif
