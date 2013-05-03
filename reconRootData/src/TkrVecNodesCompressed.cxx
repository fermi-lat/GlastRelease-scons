#include "reconRootData/TkrVecNodesCompressed.h"
#include "Riostream.h"

#include "ReconObjectManager.h"

ClassImp(TkrVecNodeCompressed)

void TkrVecNodeCompressed::Clear(Option_t* /* option */) {
    // Make sure everything is properly zeroed
    m_topPointClusterX   = 0;
    m_topPointClusterY   = 0;
    m_topPointStatus     = 0;
    m_botPointClusterX   = 0;
    m_botPointClusterY   = 0;
    m_botPointStatus     = 0;
    m_linkStatus         = 0;
    m_parent             = 0;
    m_statusBits         = 0;
    m_rmsAngleSum        = 0.;
    m_numAnglesInSum     = 0.;
    m_leaves             = 0;
    m_branches           = 0;
    m_depth              = 0;
    m_bestNumBiLayers    = 0;
    m_bestRmsAngle       = 0.;

    // Don't forget to zap the array of pointers to daughters!
    TObjArray::Clear();
}

void TkrVecNodeCompressed::Print(Option_t *option) const {
    TObject::Print(option);
}

void TkrVecNodeCompressed::initializeInfo(TkrVecNodeCompressed* parent, 
                                          UInt_t                statusBits,
                                          Double_t              rmsAngleSum,
                                          Int_t                 numAnglesInSum,
                                          Int_t                 leaves,
                                          Int_t                 branches,
                                          Int_t                 depth,
                                          Int_t                 bestNumBiLayers,
                                          Double_t              bestRmsAngle)
{
    m_parent          = parent;
    m_statusBits      = statusBits;
    m_rmsAngleSum     = rmsAngleSum;
    m_numAnglesInSum  = numAnglesInSum;
    m_leaves          = leaves;
    m_branches        = branches;
    m_depth           = depth;
    m_bestNumBiLayers = bestNumBiLayers;
    m_bestRmsAngle    = bestRmsAngle;
}
void TkrVecNodeCompressed::initializeLinkInfo(TkrCluster* topPointClusterX,
                                              TkrCluster* topPointClusterY,
                                              UInt_t      topPointStatus,
                                              TkrCluster* botPointClusterX,
                                              TkrCluster* botPointClusterY,
                                              UInt_t      botPointStatus,
                                              UInt_t      linkStatus)
{
    m_topPointClusterX = topPointClusterX;
    m_topPointClusterY = topPointClusterY;
    m_topPointStatus   = topPointStatus;
    m_botPointClusterX = botPointClusterX;
    m_botPointClusterY = botPointClusterY;
    m_botPointStatus   = botPointStatus;
    m_linkStatus       = linkStatus;
}


//======================================================
// For Unit Tests
//======================================================

void TkrVecNodeCompressed::Fake( Int_t /*ievent*/, UInt_t /*rank*/, Float_t /*randNum*/ ) {

}

#define COMPARE_IN_RANGE(att) rootdatautil::CompareInRange(get ## att(),ref.get ## att(),#att)

Bool_t TkrVecNodeCompressed::CompareInRange( const TkrVecNodeCompressed & /*ref*/, const std::string & name ) const {

    bool result = true ;

    //result = COMPARE_IN_RANGE(Energy) && result ;

    if (!result) {
        if ( name == "" ) {
            std::cout<<"Comparison ERROR for "<<ClassName()<<std::endl ;
        }
        else {
            std::cout<<"Comparison ERROR for "<<name<<std::endl ;
        }
    }
    return result ;

}

void TkrVecNodeCompressed::setTreeStartLayer(int startLayer)
{
    m_statusBits = (m_statusBits & ~START_BILAYER_BITS) 
                 | (startLayer & START_BILAYER_BITS);
    return;
}

void TkrVecNodeCompressed::setRmsAngleInfo(double rmsAngSum, int num)   
{
    m_rmsAngleSum    = rmsAngSum;
    m_numAnglesInSum = num;

    return;
}

const Double_t TkrVecNodeCompressed::getRmsAngle() const
{
    double rmsAngle = 0.;
    
    if (m_numAnglesInSum > 0) rmsAngle = sqrt(m_rmsAngleSum / m_numAnglesInSum);

    return rmsAngle;
}
    
const Double_t TkrVecNodeCompressed::getBestRmsAngle() const
{
    Double_t rmsAngle = m_bestRmsAngle; // * double(getNumBiLayers()) / double(m_depth);

    return rmsAngle;
}
    
const unsigned char TkrVecNodeCompressed::getTreeStartLayer() const
{
    unsigned char startLayer = m_statusBits & START_BILAYER_BITS;

    return startLayer;
}

const unsigned char TkrVecNodeCompressed::getCurrentBiLayer()  const
{
    unsigned char currentLayer = (m_statusBits & CURRENT_BILAYER_BITS) >> 5;

    return currentLayer;
}

const unsigned char TkrVecNodeCompressed::getBiLyrs2MainBrch()  const
{
    unsigned char currentLayer = (m_statusBits & TO_MAIN_BRANCH_BITS) >> 10;

    return currentLayer;
}

const unsigned char TkrVecNodeCompressed::getTreeId() const
{
    unsigned char treeId = (m_statusBits & TREE_ID_BITS) >> 16;

    return treeId;
}

void TkrVecNodeCompressed::setTreeId(int treeId)
{
    m_statusBits = (m_statusBits & ~TREE_ID_BITS) | ((treeId << 16) & TREE_ID_BITS);
}

void TkrVecNodeCompressed::setBiLyrs2MainBrch(int biLayersToMainBranch)
{
    m_statusBits = (m_statusBits & ~TO_MAIN_BRANCH_BITS) | ((biLayersToMainBranch << 10) & TO_MAIN_BRANCH_BITS);
}


void* TkrVecNodeCompressed::operator new(size_t /*size*/)
{
    TkrVecNodeCompressed* temp = ReconObjectManager::getPointer()->getNewTkrVecNodeCompressed();

    return temp;
}

void* TkrVecNodeCompressed::operator new(size_t /*size*/, void* vp)
{
    return vp;
}

void TkrVecNodeCompressed::operator delete(void* /*p*/)
{
    // Since we let ReconObjectManager handles memory, nothing to do here
    return;
}

