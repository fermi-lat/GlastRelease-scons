/**  @file BackgroundSelection.cxx
    @brief implementation of class BackgroundSelection
    
  $Header$  
*/

#include "BackgroundSelection.h"

#include "TFile.h"
#include "TTree.h"
#include "TLeaf.h"

#include "CLHEP/Random/RandFlat.h"

#include <stdexcept>
#include <cassert>
#include <cmath>
#include <vector>

//------------------------------------------------------------------------
BackgroundSelection::BackgroundSelection(const std::string& varname, 
                                         const std::string& rootFileDirectory,
        std::vector<std::string> disableList,
        TTree* outputTree)
: m_varname(varname)
, m_varleaf(0)
, m_event(0)
, m_eventOffset(0)
, m_outputTree(outputTree)
, m_inputTree(0)
, m_inputFile(0)
, m_disableList(disableList)
, m_rootFileDirectory(rootFileDirectory)
, m_singleFileMode(false)
{
    // find leaf corresponding to varname in the output tree
    m_varleaf = m_outputTree->GetLeaf(varname.c_str());
    if( m_varleaf==0 ) {
        throw std::runtime_error("BackgroundSelection: could not find variable "+varname);
    }
    
    setCurrentTree();

}
//------------------------------------------------------------------------
BackgroundSelection::~BackgroundSelection()
{
}

//------------------------------------------------------------------------
double BackgroundSelection::value()const
{
    return m_varleaf->GetValue();
}
//------------------------------------------------------------------------
void BackgroundSelection::selectEvent()
{

    // this is necessary due to the poor design of ROOT :-(
    TDirectory *saveDir = gDirectory;

    // make sure we have the right tree selected for new value
    setCurrentTree();    

    // check event offset for overflow:
    if (m_eventOffset >= m_inputTree->GetEntries())
        m_eventOffset = 0;

    // grab the event:
    m_inputTree->GetEvent(m_eventOffset++);
    // this is necessary due to the poor design of ROOT :-(
    saveDir->cd();
}

//------------------------------------------------------------------------
double BackgroundSelection::triggerRate()const
{
    double x(value());

    ///@todo: replace with lookup
    return 1000;
}

//------------------------------------------------------------------------
double BackgroundSelection::downlinkRate( )const
{
    double  x(value());

    ///@todo: replace with lookup
    return 100;
}
//------------------------------------------------------------------------
void BackgroundSelection::disableBranches(TTree* t) 
{
    std::vector<std::string>::const_iterator it= m_disableList.begin();
    for( ; it != m_disableList.end();  ++it) {
        t->SetBranchStatus((*it).c_str(), false);
    }

    // restore this, if needed, for cuts
    t->SetBranchStatus("FT1ZenithTheta", 1);
    t->SetBranchStatus("McSourceId", 1);
}
//------------------------------------------------------------------------
void BackgroundSelection::setCurrentTree() 
{
    if( m_singleFileMode && m_inputTree!=0 ) return; // all set up

    std::string file_name("test.root");
    std::string tree_name("MeritTuple");

    if( ! m_singleFileMode && m_rootFileDirectory.find(".root")!=std::string::npos) { 

        // if is is root file, will set it up only once
        m_singleFileMode = true;
        file_name = m_rootFileDirectory;

    }else {

        double x(value());

        ///@todo: set up lookup table for rates, files to sample from, based on currrent value
        // set file_name to be opened (and maybe close previous one?)
    }

    // open file for reading:
    m_inputFile = new TFile(file_name.c_str(), "readonly");
    if (!m_inputFile->IsOpen() ) {
        TString error = "Did not find file[" + file_name + "]";
        throw std::invalid_argument(error.Data());
    }

    // get the tree:
    m_inputTree =  dynamic_cast<TTree*>(m_inputFile->Get(tree_name.c_str()));
    if (0 == m_inputTree) {
        TString error = "Did not find tree[" + tree_name + "] in root file";
        throw std::invalid_argument(error.Data());
    }

    // start at a random location in the tree:
    m_eventOffset = (unsigned int)(RandFlat::shoot()*(m_inputTree->GetEntries() - 1));

    // point tree to buffer for copying events:
    setLeafPointers(m_inputTree);
}
//------------------------------------------------------------------------
void BackgroundSelection::setLeafPointers(TTree* pTree)
{

    // first disable the branches
    disableBranches(pTree);

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
            std::stringstream msg;
            msg << "BackgroundSelection::setLeafPointesrs: did not find the leaf " << leaf->GetName();
            std::cerr << msg.str() << std::endl;
            throw std::runtime_error( msg.str());
        }else{
            leaf->SetAddress(otherleaf->GetValuePointer());
        }
    }
}

