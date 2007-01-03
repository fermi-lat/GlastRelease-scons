/**  @file BackgroundSelection.cxx
    @brief implementation of class BackgroundSelection
    
  $Header$  
*/

#include "BackgroundSelection.h"
#include "XmlFetchEvents.h"

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TLeaf.h"

#include "CLHEP/Random/RandFlat.h"

#include <stdexcept>
#include <cassert>
#include <cmath>
#include <vector>

//------------------------------------------------------------------------
BackgroundSelection::BackgroundSelection(const std::string& varname, 
                                         const std::string& xmlFileName,
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
, m_fetch(0)
, m_observer(0)
{
    // find leaf corresponding to varname in the output tree
    m_varleaf = m_outputTree->GetLeaf(varname.c_str());
    if( m_varleaf==0 ) {
        throw std::runtime_error("BackgroundSelection: could not find variable "+varname);
    }

    // If passed an XML file, set up the XmlFetchEvents object
    if (xmlFileName.find(".xml")!=std::string::npos){
        m_fetch = new XmlFetchEvents(xmlFileName, varname);
    }else{ // did not find an option
        throw std::invalid_argument("BackgroundSelection: expect to find an XML file");
    }

    //---------------------------------------------------
    // ----------- for getting notification -------------
    class Observer : public TObject {
    public:
        Observer(BackgroundSelection& bg): m_bg(bg){}
        virtual bool Notify(){m_bg.notify(); return true;} // virtual in TObject
    private:
        BackgroundSelection& m_bg;
    };
    m_observer = new Observer(*this);

    // prime the pump with the mean of the limits for initial downlink, trigger rates
    double mean( (m_fetch->minVal() + m_fetch->maxVal())/2.);

    setCurrentTree(mean);
}
//------------------------------------------------------------------------
BackgroundSelection::~BackgroundSelection()
{
    delete m_fetch;
    delete m_observer;
}
//------------------------------------------------------------------------
void BackgroundSelection::notify()
{
    std::cout << " Got notified " << std::endl;
}

//------------------------------------------------------------------------
double BackgroundSelection::value()const
{
    double x(m_varleaf->GetValue());
    return x;
}
//------------------------------------------------------------------------
void BackgroundSelection::selectEvent()
{


    // make sure we have the right tree selected for new value
    setCurrentTree(value());    

    // grab the next event, cycling if needed
    Long64_t bytes = m_inputTree->LoadTree(m_eventOffset++);
    if( bytes<=0){ m_eventOffset = 0; 
        bytes = m_inputTree->LoadTree(m_eventOffset++);
    }
    if( bytes<=0) throw std::runtime_error("BackgroundSelection::selectEvent -- could not read file");

}
//------------------------------------------------------------------------
void BackgroundSelection::disableBranches(TTree* t) 
{
    std::vector<std::string>::const_iterator it= m_disableList.begin();
    for( ; it != m_disableList.end();  ++it) {
        t->SetBranchStatus((*it).c_str(), false);
    }

}
//------------------------------------------------------------------------
void BackgroundSelection::setCurrentTree(double x) 
{
    TTree * nextTree(0);
    // this is necessary due to the design of ROOT :-(
    TDirectory *saveDir = gDirectory;

    TChain* chain = new TChain();
    int stat= m_fetch->getFiles(x, chain);
    if( stat!=0 ) throw std::runtime_error("BackgaroundSelection::setCurrentTree: invalid tree");
    chain->SetNotify(m_observer);
    m_inputTree = chain;
    // this is necessary due to the design of ROOT :-(
    saveDir->cd();


    if( nextTree==0) return; // still using the same treee

    // a new tree: must set it up for sequential access from random point
    m_inputTree = nextTree;

    // start at a random location in the tree:
    double length (m_inputTree->GetEntries());
    if( length==0 ) {
        throw std::runtime_error("BackgroundSelection::setCurrentTree: no events in the tree");
    }
    
    m_eventOffset = (unsigned int)(RandFlat::shoot()*(length - 1));

    // point tree to buffer for copying events:
    setLeafPointers(m_inputTree);
    m_triggerRate = m_fetch->getAttributeValue("triggerRate", x);
    m_downlinkRate=m_fetch->getAttributeValue("downlinkRate", x);

    int ret =m_inputTree->GetEntry(0); // make sure this works
    if( ret<0 ) throw std::runtime_error("BackgroundSelection::setCurrentTree: could not read");

}
//------------------------------------------------------------------------
void BackgroundSelection::setLeafPointers(TTree* pTree)
{
    // first disable the branches
    disableBranches(pTree);

    // iteration over active leaves -- copied from Ttree::Show()

    TObjArray *leaves  = pTree->GetListOfLeaves();
    if( leaves==0 ) throw std::runtime_error("BackgroundSelection::setLeafPointers: invalid tree, no leaves found");

    Int_t nleaves = leaves->GetEntries();
    if (nleaves==0) throw std::runtime_error("BackgroundSelection::setLeafPointers: invalid tree, no leaves found");
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
            //std::cout << "setting pointer for leaf " << leaf->GetName() << std::endl;
            leaf->SetAddress(otherleaf->GetValuePointer());
        }
    }
}

