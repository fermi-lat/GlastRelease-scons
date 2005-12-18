/**  @file BackgroundSelection.cxx
    @brief implementation of class BackgroundSelection
    
    
*/

#include "BackgroundSelection.h"

#include "TFile.h"
#include "TTree.h"

#include <stdexcept>

//------------------------------------------------------------------------
BackgroundSelection::BackgroundSelection(const std::string& filename, const std::string& treename)
: m_file(0)
, m_event(0)
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
    // TODO: fill this
    // note that we do not want to copy any leaf value if the corresponding branch is inactive

}

//------------------------------------------------------------------------
TTree* BackgroundSelection::selectEvent(double /* maglat */)
{
    // TODO: have this depend on magnetic latitude
    if( ++m_event > m_tree->GetEntries() ) {
        m_event = 0;
    }
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
