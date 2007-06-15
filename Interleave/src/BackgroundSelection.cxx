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
#include <sstream>
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
    double mean( (m_fetch->minValFullRange() + m_fetch->maxValFullRange())/2.);

    setCurrentTree(mean);
}
//------------------------------------------------------------------------
BackgroundSelection::~BackgroundSelection()
{
    delete m_fetch;
    delete m_observer;
    delete m_inputTree;
}
//------------------------------------------------------------------------
void BackgroundSelection::notify()
{
    setLeafPointers();
}

const std::string& BackgroundSelection::sourceName() const 
{
    return m_fetch->name();
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

    double x(value());
#if 1 // have to do this since my specification was not understood
    x = fabs(x); 
#endif
    if( !m_fetch->isValid(x) ){
        std::stringstream msg;
        msg << "BackgaroundSelection::selectEvent: called with "<< name() 
            <<" = "<< x << " is not in range " 
            <<  m_fetch->minValFullRange() << ", to " 
            <<  m_fetch->maxValFullRange() << std::endl;
        throw std::runtime_error(msg.str());
    }

    // make sure we have the right tree selected for new value    // if still valid, do not change
    if( !m_fetch->isCurrent(x) ){
        setCurrentTree(x);   
    }

    // grab the next event, cycling if needed
    //Long64_t code = m_inputTree->LoadTree(m_eventOffset++);
    int code = m_inputTree->GetEntry(m_eventOffset++);
    if( code<0){ m_eventOffset = 0; 
        code = m_inputTree->LoadTree(m_eventOffset++);
        if( code < 0 ){
            throw std::runtime_error("BackgroundSelection::selectEvent -- could not read file");
        }
    }
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

    // replace the TChain
    delete m_inputTree; 
    m_inputTree = new TChain("MeritTuple");
    //m_inputTree->SetNotify(m_observer);

    // this is necessary due to the design of ROOT :-(
    TDirectory *saveDir = gDirectory;
    
    int stat= m_fetch->getFiles(x, dynamic_cast<TChain*>(m_inputTree));
    if( stat!=0 ){
        std::stringstream msg;
        msg << "BackgaroundSelection::setCurrentTree: called with "<< name() 
            <<" = "<< x << ", getFiles returned error code " << stat;
        throw std::runtime_error(msg.str());
    }

    // this is necessary due to the design of ROOT :-(
    saveDir->cd();

    // start at a random location in the tree:
    double length (m_inputTree->GetEntries());
    if( length==0 ) {
        throw std::runtime_error("BackgroundSelection::setCurrentTree: no events in the tree");
    }
    
    m_eventOffset = (unsigned int)(RandFlat::shoot()*(length - 1));

    // point tree to buffer for copying events:
    setLeafPointers();
    m_triggerRate = m_fetch->getAttributeValue("triggerRate", x);
    m_downlinkRate=m_fetch->getAttributeValue("downlinkRate", x);

    //int ret =m_inputTree->LoadTree(0); // make sure this works
    int ret = m_inputTree->GetEntry(0);
    if( ret<0 ) throw std::runtime_error("BackgroundSelection::setCurrentTree: could not read");

}
//------------------------------------------------------------------------
void BackgroundSelection::setLeafPointers()
{
    // first disable the branches
    disableBranches(m_inputTree);

    // iteration over active leaves -- copied from Ttree::Show()

    TObjArray *leaves  = m_inputTree->GetListOfLeaves();
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
            msg << "BackgroundSelection::setLeafPointers: did not find the leaf " << leaf->GetName();
            std::cerr << msg.str() << std::endl;
            throw std::runtime_error( msg.str());
        }else{
            //std::cout << "setting pointer for leaf " << leaf->GetName() << std::endl;
            leaf->SetAddress(otherleaf->GetValuePointer());
        }
    }
}

