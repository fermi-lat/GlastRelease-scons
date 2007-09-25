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
    TFile* check = nextSkim->openFile(aFile->GetName());
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
  Long64_t* eventIndex = new Long64_t(-1);
  Long64_t* fileSetIndexOut  = new Long64_t(-1);
  Long64_t* fileSetOffset = new Long64_t(0);
  TTree* eventTreeOut = new TTree("CompositeEvents","CompositeEvents");
  eventTreeOut->Branch("Event_Index",(void*)eventIndex,"Event_Index/L");
  eventTreeOut->Branch("Event_FileSetIndex",(void*)fileSetIndexOut,"Event_FileSetIndex/L");
  eventTreeOut->Branch("Event_FileSetOffset",(void*)fileSetOffset,"Event_FileSetOffset/L");

  // Loop on cels
  TIterator* itr = cels.MakeIterator();
  TObject* aSkimObj(0);
  while ( (aSkimObj = itr->Next()) != 0 ) {
    CompositeEventList* aSkim = dynamic_cast<CompositeEventList*>(aSkimObj);
    if ( 0 == aSkim ) {
      std::cerr << "Input object not a CompositeEventList " << aSkimObj->GetName() << std::endl;
      return 0;
    }    
    
    // Add the event & file trees to their respective lists
    entryTreeList.AddLast(aSkim->entryTree());
    fileTreeList.AddLast(aSkim->fileTree());    

    Long64_t fileSetIndexIn(0);
    Long64_t fileSetIndexSave(-1);
    TTree* eventTree = aSkim->linkTree();
    eventTree->SetBranchAddress("Event_FileSetIndex",&fileSetIndexIn);

    // Redo the data in the Event tree
    Long64_t nEvtCurrent = eventTree->GetEntries();
    for ( Long64_t iEvt(0); iEvt < nEvtCurrent; iEvt++ ) {
      Int_t nB = eventTree->GetEntry(iEvt);
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
	delete eventTree;
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
      eventTreeOut->Fill();
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
  CompositeEventList * theCel = new CompositeEventList(*eventTreeOut,*fileTree,*entryTree) ;

  // Write the Merged File, if requested
  if ( outFile != 0 )
   {
    outFile->Write();
    outFile->Close();
   }
  return theCel ;  
 }


#endif
