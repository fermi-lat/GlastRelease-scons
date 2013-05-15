/** @file TreeClusterRelation.h

* @class TreeClusterRelation
*
* @brief This is the root version of the TreeClusterRelation object. The TreeClusterRelation object
*        serves to relate a Tree object in the Tracker to a Cluster object in the Calorimeter
*
* last modified 02/01/2013
*
* @authors Usher
*
* $Header: 
*/

#ifndef ROOT__TreeClusterRelation_H
#define ROOT__TreeClusterRelation_H

#include "TObject.h"
#include "TRef.h"
#include "reconRootData/TkrTree.h"
#include "reconRootData/CalCluster.h"

class TreeClusterRelation: public TObject
{
public:
    // constructors - null case
    TreeClusterRelation() {Clear();}

    // constructurs - have the information case
//    TreeClusterRelation(TreeClusterRelation* parent, const TkrVecPointsLink* associatedLink);

    // destructor
    virtual ~TreeClusterRelation() {Clear();}

    // Root mandated methods
    void Clear(Option_t *option="");
    void Fake( Int_t ievent, UInt_t rank, Float_t randNum ) ; // for tests
    Bool_t CompareInRange( const TreeClusterRelation &, const std::string & name = "" ) const ; // for tests
    void Print(Option_t *option="") const;

    // Initializer
    void initializeInfo(TkrTree*    tree,
                        CalCluster* cluster,
                        Double_t    treeClusDoca,
                        Double_t    treeClusCosAngle,
                        Double_t    treeClusDistAtZ,
                        Double_t    clusEnergy);

    /// Status word bits organized like:
    ///        |  0   0   0   0  |  0   0   0   0  |  0   0   0   0  |  0   0   0   0   |
    ///         <               > <               > <               >  <              >

    /// Set data members
    void                     setTree(TkrTree* tree)                         {m_tree             = tree;}
    void                     setCluster(CalCluster* cluster)                {m_cluster          = cluster;}
    void                     setTreeClusDoca(Double_t treeClusDoca)         {m_treeClusDoca     = treeClusDoca;}
    void                     setTreeClusCosAngle(Double_t treeClusCosAngle) {m_treeClusCosAngle = treeClusCosAngle;}
    void                     setTreeClusDistAtZ(Double_t treeClusDistAtZ)   {m_treeClusDistAtZ  = treeClusDistAtZ;}
    void                     setClusEnergy(Double_t energy)                 {m_clusEnergy       = energy;}

    /// Recover data members
    Double_t                 getTreeClusDoca()                              {return m_treeClusDoca;}
    Double_t                 getTreeClusCosAngle()                          {return m_treeClusCosAngle;}
    Double_t                 getTreeClusDistAtZ()                           {return m_treeClusDistAtZ;}
    Double_t                 getClusEnergy()                                {return m_clusEnergy;}

    const TkrTree*           getTree()                                const {return (const TkrTree*)m_tree.GetObject();}
    const CalCluster*        getCluster()                             const {return (const CalCluster*)m_cluster.GetObject();}
    const Double_t           getTreeClusDoca()                        const {return m_treeClusDoca;}
    const Double_t           getTreeClusCosAngle()                    const {return m_treeClusCosAngle;}
    const Double_t           getTreeClusDistAtZ()                     const {return m_treeClusDistAtZ;}
    const Double_t           getClusEnergy()                          const {return m_clusEnergy;}

private:
    TRef       m_tree;
    TRef       m_cluster;
    Double32_t m_treeClusDoca;
    Double32_t m_treeClusCosAngle;
    Double32_t m_treeClusDistAtZ;
    Double32_t m_clusEnergy;

    ClassDef(TreeClusterRelation,2)
};


#endif // __TreeClusterRelation_H
