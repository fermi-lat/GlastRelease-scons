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

//------------------------------------------------------------------------
BackgroundSelection::BackgroundSelection(const std::string& filename, const std::string& treename)
: m_event(0)
{

    TFile* m_file = new TFile(filename.c_str(), "readonly");
    if( 0==m_file || !m_file->IsOpen()){
        throw std::invalid_argument( "Could not open the root file ");
    }

    m_tree =  dynamic_cast<TTree*>(m_file->Get(treename.c_str()));
    if( 0==m_tree) {
        throw std::invalid_argument("Did not find the tree in the root file " );
    }
}
//------------------------------------------------------------------------
BackgroundSelection::~BackgroundSelection()
{
    delete m_file;
}
//------------------------------------------------------------------------
void BackgroundSelection::copyTreeData( TTree* out)
{

}

//------------------------------------------------------------------------
void BackgroundSelection::setLeafPointers(TTree* other)
{

    // iteration over active leaves -- copied from Ttree::Show()

    TObjArray *leaves  = m_tree->GetListOfLeaves();
    Int_t nleaves = leaves->GetEntriesFast();
    for (Int_t i=0;i<nleaves;i++) {
        TLeaf *leaf = (TLeaf*)leaves->UncheckedAt(i);
        TBranch *branch = leaf->GetBranch();
        if (branch->TestBit(kDoNotProcess)) continue;
        Int_t len = leaf->GetLen();
        if (len <= 0) continue;

        // here copy the value pointer from the merit tuple
        TLeaf* otherleaf = other->GetLeaf(leaf->GetName());
        if( otherleaf==0){
            std::cout << "No leaf " << leaf->GetName() << std::endl;
        }else{
            leaf->SetAddress(otherleaf->GetValuePointer());
        }
    }
}

//------------------------------------------------------------------------
TTree* BackgroundSelection::selectEvent(double /* maglat */)
{
    // TODO: have this depend on magnetic latitude
    if( ++m_event > m_tree->GetEntries() ) {
        m_event = 0;
    }
    // note: if this overrites the other TTree, a GetEvent should be deferred
    m_tree->GetEvent(m_event);
    return m_tree;
}

//------------------------------------------------------------------------
double BackgroundSelection::triggerRate(double/* maglat */)
{
    //TODO
    return 20.; // testing
}

//------------------------------------------------------------------------
double BackgroundSelection::downlinkRate(double/* maglat */)
{
    //TODO
    return 1.; //testing
}
void BackgroundSelection::disable(const char* pattern)
{
    m_tree->SetBranchStatus( pattern, 0);
}
