/**  @file BackgroundSelection.cxx
    @brief implementation of class BackgroundSelection
    
  $Header$  
*/

#include "BackgroundSelection.h"

#include "TFile.h"
#include "TTree.h"
#include "TLeaf.h"

#include <stdexcept>
#include <cassert>
#include <math.h>

//------------------------------------------------------------------------
BackgroundSelection::BackgroundSelection(const std::string& filename, 
                                         const std::string& treename
                                         , TTree* outputTree)
: m_event(0)
, m_tree(0)
, m_file(0)
, m_outputTree(outputTree)
{

    m_file = new TFile(filename.c_str(), "readonly");
    if( 0==m_file || !m_file->IsOpen()){
        throw std::invalid_argument( "Could not open the root file ");
    }

    for (int latBin=-42; latBin<42; latBin++) {
      int binIndex = latBin + 42;
      m_inputTreeIndexes[binIndex] = 0;
      TString tree_name;
      tree_name += latBin;
      tree_name += "to";
      tree_name += (latBin + 1);
      m_inputTrees[binIndex] = dynamic_cast<TTree*>(m_file->Get(tree_name.Data()));
      if (0 == m_inputTrees[binIndex]) {
	TString error = "Did not find tree[" + tree_name + "] in root file";
        throw std::invalid_argument(error.Data());
      }
    }
}
//------------------------------------------------------------------------
BackgroundSelection::~BackgroundSelection()
{
    for (int latBin=-42; latBin<42; latBin++) {
      int binIndex = latBin + 42;
      if (m_inputTrees[binIndex] != 0)
	delete m_inputTrees[binIndex];
    }
    delete m_file;
}

//------------------------------------------------------------------------
void BackgroundSelection::setLeafPointers(TTree* pTree)
{

    // iteration over active leaves -- copied from Ttree::Show()

    TObjArray *leaves  = pTree->GetListOfLeaves();
    Int_t nleaves = leaves->GetEntriesFast();
    for (Int_t i=0;i<nleaves;i++) {
        TLeaf *leaf = (TLeaf*)leaves->UncheckedAt(i);
        TBranch *branch = leaf->GetBranch();
        if (branch->TestBit(kDoNotProcess)) continue;
        Int_t len = leaf->GetLen();
        if (len <= 0) continue;

        // here copy the value pointer from the merit tuple
        TLeaf* otherleaf = m_outputTree->GetLeaf(leaf->GetName());
        if( otherleaf==0){
            std::cout << "No leaf " << leaf->GetName() << std::endl;
        }else{
            leaf->SetAddress(otherleaf->GetValuePointer());
        }
    }
}

//------------------------------------------------------------------------
void BackgroundSelection::selectEvent(double maglat )
{
    // TODO: have this depend on magnetic latitude

    int binIndex = (int)(floor(maglat)) + 42;

    TTree* pTree = m_inputTrees[binIndex];
    if (m_inputTreeIndexes[binIndex] >= pTree->GetEntries())
      m_inputTreeIndexes[binIndex] = 0;

    Long64_t nEvent = m_inputTreeIndexes[binIndex++];
    
    setLeafPointers(pTree); // may only have to be done once

    pTree->GetEvent(nEvent);
}

//------------------------------------------------------------------------
double BackgroundSelection::triggerRate(double/* maglat */)
{
    //TODO
    return 2000.; // testing
}

//------------------------------------------------------------------------
double BackgroundSelection::downlinkRate(double/* maglat */)
{
    //TODO
    return 300.; //testing
}
void BackgroundSelection::disable(const char* pattern)
{
    m_tree->SetBranchStatus( pattern, 0);
}
