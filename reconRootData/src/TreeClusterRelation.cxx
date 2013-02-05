#include "reconRootData/TreeClusterRelation.h"
#include "Riostream.h"

#include "ReconObjectManager.h"

ClassImp(TreeClusterRelation)
void TreeClusterRelation::initializeInfo(TkrTree*    tree,
                                         CalCluster* cluster,
                                         Double_t    treeClusDoca,
                                         Double_t    treeClusCosAngle,
                                         Double_t    treeClusDistAtZ,
                                         Double_t    clusEnergy)
{
    m_tree             = tree;
    m_cluster          = cluster;
    m_treeClusDoca     = treeClusDoca;
    m_treeClusCosAngle = treeClusCosAngle;
    m_treeClusDistAtZ  = treeClusDistAtZ;
    m_clusEnergy       = clusEnergy;
}
    
void TreeClusterRelation::Clear(Option_t* /* option */) 
{
    // Make sure everything is properly zeroed
    m_tree             = 0;
    m_cluster          = 0;
    m_treeClusDoca     = 0.;
    m_treeClusCosAngle = 0.;
    m_treeClusDistAtZ  = 0.;
    m_clusEnergy       = 0.;
}

void TreeClusterRelation::Print(Option_t *option) const {
    TObject::Print(option);
}



//======================================================
// For Unit Tests
//======================================================

void TreeClusterRelation::Fake( Int_t /*ievent*/, UInt_t /*rank*/, Float_t /*randNum*/ ) {

}

#define COMPARE_IN_RANGE(att) rootdatautil::CompareInRange(get ## att(),ref.get ## att(),#att)

Bool_t TreeClusterRelation::CompareInRange( const TreeClusterRelation & /*ref*/, const std::string & name ) const {

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
