#include "reconRootData/TkrVecNodes.h"
#include "Riostream.h"

#include "ReconObjectManager.h"

ClassImp(TkrVecNode)

void TkrVecNode::Clear(Option_t* /* option */) {
    // Make sure everything is properly zeroed
    m_associatedLink  = 0;
    m_parent          = 0;
    m_statusBits      = 0;
    m_rmsAngleSum     = 0.;
    m_numAnglesInSum  = 0.;
    m_leaves          = 0;
    m_branches        = 0;
    m_depth           = 0;
    m_bestNumBiLayers = 0;
    m_bestRmsAngle    = 0.;

    // Don't forget to zap the array of pointers to daughters!
    TRefArray::Clear();
}

void TkrVecNode::Print(Option_t *option) const {
    TObject::Print(option);
}

void TkrVecNode::initializeInfo(TkrVecNode*             parent, 
                                const TkrVecPointsLink* associatedLink, 
                                UInt_t                  statusBits,
                                Double_t                rmsAngleSum,
                                Int_t                   numAnglesInSum,
                                Int_t                   leaves,
                                Int_t                   branches,
                                Int_t                   depth,
                                Int_t                   bestNumBiLayers,
                                Double_t                bestRmsAngle)
{
    m_associatedLink  = associatedLink;
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


//======================================================
// For Unit Tests
//======================================================

void TkrVecNode::Fake( Int_t /*ievent*/, UInt_t /*rank*/, Float_t /*randNum*/ ) {

}

#define COMPARE_IN_RANGE(att) rootdatautil::CompareInRange(get ## att(),ref.get ## att(),#att)

Bool_t TkrVecNode::CompareInRange( const TkrVecNode & /*ref*/, const std::string & name ) const {

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

void TkrVecNode::setTreeStartLayer(int startLayer)
{
    m_statusBits = (m_statusBits & ~START_BILAYER_BITS) 
                 | (startLayer & START_BILAYER_BITS);
    return;
}

void TkrVecNode::setRmsAngleInfo(double rmsAngSum, int num)   
{
    m_rmsAngleSum    = rmsAngSum;
    m_numAnglesInSum = num;

    return;
}

const Double_t TkrVecNode::getRmsAngle() const
{
    double rmsAngle = 0.;
    
    if (m_numAnglesInSum > 0) rmsAngle = sqrt(m_rmsAngleSum / m_numAnglesInSum);

    return rmsAngle;
}
    
const Double_t TkrVecNode::getBestRmsAngle() const
{
    Double_t rmsAngle = m_bestRmsAngle; // * double(getNumBiLayers()) / double(m_depth);

    return rmsAngle;
}
    
const unsigned char TkrVecNode::getTreeStartLayer() const
{
    unsigned char startLayer = m_statusBits & START_BILAYER_BITS;

    return startLayer;
}

const unsigned char TkrVecNode::getCurrentBiLayer()  const
{
    unsigned char currentLayer = (m_statusBits & CURRENT_BILAYER_BITS) >> 5;

    return currentLayer;
}

const unsigned char TkrVecNode::getBiLyrs2MainBrch()  const
{
    unsigned char currentLayer = (m_statusBits & TO_MAIN_BRANCH_BITS) >> 10;

    return currentLayer;
}

const unsigned char TkrVecNode::getTreeId() const
{
    unsigned char treeId = (m_statusBits & TREE_ID_BITS) >> 16;

    return treeId;
}

void TkrVecNode::setTreeId(int treeId)
{
    m_statusBits = (m_statusBits & ~TREE_ID_BITS) | ((treeId << 16) & TREE_ID_BITS);
}

void TkrVecNode::setBiLyrs2MainBrch(int biLayersToMainBranch)
{
    m_statusBits = (m_statusBits & ~TO_MAIN_BRANCH_BITS) | ((biLayersToMainBranch << 10) & TO_MAIN_BRANCH_BITS);
}


void* TkrVecNode::operator new(size_t /*size*/)
{
    TkrVecNode* temp = ReconObjectManager::getPointer()->getNewTkrVecNode();

    return temp;
}

void* TkrVecNode::operator new(size_t /*size*/, void* vp)
{
    return vp;
}

void TkrVecNode::operator delete(void* /*p*/)
{
    // Since we let ReconObjectManager handles memory, nothing to do here
    return;
}

