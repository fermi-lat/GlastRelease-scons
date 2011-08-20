/** @file TreeCalClusterAssociator.h
 * @class TreeCalClusterAssociator
 *
 * @brief Provides X-Y space points from the same tower from TkrClusterCol
 *
 * @author Tracy Usher
 *
 * $Header$
 *
*/

#ifndef __TreeClusterRelation_H
#define __TreeClusterRelation_H 1

#include "GaudiKernel/ObjectVector.h"
#include "GaudiKernel/ContainedObject.h"

//#include <vector>

#include "Event/Recon/TkrRecon/TkrTree.h"
#include "Event/Recon/CalRecon/CalCluster.h"

// Declare Gaudi object interface ID
static const CLID& CLID_TreeClusterRelation = InterfaceID("TreeClusterRelation",  1, 0);

// Enclose the class definition in the Event namespace
namespace Event
{
    class TreeClusterRelation: virtual public ContainedObject
    {
    public:
        TreeClusterRelation() : m_tree(0),
                                m_cluster(0),
                                m_treeClusDoca(0.),
                                m_treeClusCosAngle(0.),
                                m_treeClusDistAtZ(0.),
                                m_clusEnergy(0.)
        {}

        TreeClusterRelation(Event::TkrTree* tree,
                            Event::CalCluster* cluster,
                            double             treeClusDoca,
                            double             treeClusCosAngle,
                            double             treeClusDistAtZ,
                            double             clusEnergy)
                            : m_tree(tree),
                              m_cluster(cluster),
                              m_treeClusDoca(treeClusDoca),
                              m_treeClusCosAngle(treeClusCosAngle),
                              m_treeClusDistAtZ(treeClusDistAtZ),
                              m_clusEnergy(clusEnergy)
        {}

       ~TreeClusterRelation() {}

        //! Retrieve pointer to class defininition structure
        virtual const CLID& clID()    const   { return TreeClusterRelation::classID(); }
        static  const CLID& classID()         { return CLID_TreeClusterRelation; }

       void                     setTree(Event::TkrTree* tree)                {m_tree             = tree;}
       void                     setCluster(Event::CalCluster* cluster)       {m_cluster          = cluster;}
       void                     setTreeClusDoca(double treeClusDoca)         {m_treeClusDoca     = treeClusDoca;}
       void                     setTreeClusCosAngle(double treeClusCosAngle) {m_treeClusCosAngle = treeClusCosAngle;}
       void                     setTreeClusDistAtZ(double treeClusDistAtZ)   {m_treeClusDistAtZ  = treeClusDistAtZ;}
       void                     setClusEnergy(double energy)                 {m_clusEnergy       = energy;}

       Event::TkrTree*          getTree()                                    {return m_tree;}
       Event::CalCluster*       getCluster()                                 {return m_cluster;}
       double                   getTreeClusDoca()                            {return m_treeClusDoca;}
       double                   getTreeClusCosAngle()                        {return m_treeClusCosAngle;}
       double                   getTreeClusDistAtZ()                         {return m_treeClusDistAtZ;}
       double                   getClusEnergy()                              {return m_clusEnergy;}

       const Event::TkrTree*    getTree()                              const {return m_tree;}
       const Event::CalCluster* getCluster()                           const {return m_cluster;}
       const double             getTreeClusDoca()                      const {return m_treeClusDoca;}
       const double             getTreeClusCosAngle()                  const {return m_treeClusCosAngle;}
       const double             getTreeClusDistAtZ()                   const {return m_treeClusDistAtZ;}
       const double             getClusEnergy()                        const {return m_clusEnergy;}

       const bool operator<(const TreeClusterRelation* right) const;

    private:
        Event::TkrTree*    m_tree;
        Event::CalCluster* m_cluster;
        double             m_treeClusDoca;
        double             m_treeClusCosAngle;
        double             m_treeClusDistAtZ;
        double             m_clusEnergy;
    };

    inline const bool TreeClusterRelation::operator<(const TreeClusterRelation* right) const
    {
        return true;
    }

    // Typedef the object that will be the container/owner in the TDS
    typedef ObjectVector<TreeClusterRelation>                    TreeClusterRelationCol;

    typedef std::vector<TreeClusterRelation*>                    TreeClusterRelationVec;

    typedef std::map<Event::TkrTree*,    TreeClusterRelationVec> TreeToRelationMapDef;
    typedef std::map<Event::CalCluster*, TreeClusterRelationVec> ClusterToRelationMapDef;

    // Define the TDS version of the map between trees and the relations
    static const CLID& CLID_TreeToRelationMap  = InterfaceID("TreeToRelationMap",  1, 0);

    class TreeToRelationMap : public TreeToRelationMapDef, public DataObject 
    {
        //! Retrieve pointer to class defininition structure
        virtual const CLID& clID() const   { return TreeToRelationMap::classID(); }
        static const CLID& classID()       { return CLID_TreeToRelationMap; }
    };

    // Define the TDS version of the map between trees and the relations
    static const CLID& CLID_ClusterToRelationMap  = InterfaceID("ClusterToRelationMap",  1, 0);

    class ClusterToRelationMap : public ClusterToRelationMapDef, public DataObject 
    {
        //! Retrieve pointer to class defininition structure
        virtual const CLID& clID() const   { return ClusterToRelationMap::classID(); }
        static const CLID& classID()       { return CLID_ClusterToRelationMap; }
    };
};

#endif
