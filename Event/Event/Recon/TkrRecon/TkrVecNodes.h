/**
 * @class TkrVecNodes
 *
 * @brief This defines a class to associate TkrVecPoints into "links"
 *
 * @author Tracy Usher
 *
 * $Header$
 */

#ifndef TkrVecNode_h
#define TkrVecNode_h

#include "Event/Recon/TkrRecon/TkrVecPointsLink.h"

#include "GaudiKernel/ObjectVector.h"
#include "GaudiKernel/ContainedObject.h"
#include "GaudiKernel/IInterface.h"

#include "Event/RelTable/RelTable.h"

#include <set>

// Declare Gaudi object interface ID
static const CLID& CLID_TkrVecNode = InterfaceID("TkrVecNode",  1, 0);

namespace Event {  // NameSpace

// Forward declaration
class TkrVecNode;

// Define a comparator which will be used to set the order of nodes inserted into
// the set of daughter nodes. Idea is "best" branch will be first. 
class TkrVecNodesComparator
{
public:
    // Define operator to facilitate sorting
    const bool operator()(const TkrVecNode* left, const TkrVecNode* right) const;
};

// Typedef's to define the set which will contain our daughters
//typedef std::set<TkrVecNode*, TkrVecNodesComparator> TkrVecNodeSet;
typedef std::list<TkrVecNode*> TkrVecNodeSet;

class TkrVecNode: public TkrVecNodeSet, virtual public ContainedObject
{
public:
    // enumerate the status bits
    enum StatusBits {NODE_IS_SACRED           = 0x80000000,
                     NODE_CAN_BE_SHARED       = 0x01000000,
                     NODE_ON_BEST_BRANCH      = 0x02000000,
                     NODE_ON_NEXT_BEST_BRANCH = 0x04000000};

    enum LayerMask  {START_BILAYER_BITS       = 0x0000001F,
                     CURRENT_BILAYER_BITS     = 0x000003E0,
                     TO_MAIN_BRANCH_BITS      = 0x00007C00,
                     TREE_ID_BITS             = 0x00FF0000};
    // Constructors
    TkrVecNode(TkrVecNode*             parent, 
               const TkrVecPointsLink* associatedLink);

    virtual ~TkrVecNode(); 

    //! Retrieve pointer to class defininition structure
    virtual const CLID& clID()    const   { return TkrVecNode::classID(); }
    static  const CLID& classID()         { return CLID_TkrVecNode; }

    // Methods to manipulate the data members
    // Change the parent for this node
    void changeParentNode(TkrVecNode* newParent) {m_parent          = newParent;}
    // Method to update the rms angle going TO this node
    void setRmsAngleInfo(double rms, int num);
    // Set the number of leaves below this node
    void setNumLeaves(int numLeaves)             {m_leaves          = numLeaves;}
    // Set the number of branches below this nodde
    void setNumBranches(int numBranches)         {m_branches        = numBranches;}
    // Set the depth below this node
    void setDepth(int depth)                     {m_depth           = depth;}
    // Set the number of bilayers below this node for best branch
    void setBestNumBiLayers(int numBiLayers)     {m_bestNumBiLayers = numBiLayers;}
    // Set the RMS deviation below this node for best branch
    void setBestRmsAngle(double rmsAngle)        {m_bestRmsAngle    = rmsAngle;}
    // Set the Tree ID #
    void setTreeId(int treeId);
    // Set status bits
    void setNodeSacred()                         {m_statusBits |= NODE_IS_SACRED;}
    void setNodeCanBeShared()                    {m_statusBits |= NODE_CAN_BE_SHARED;}
    void setNodeOnBestBranch()                   {m_statusBits |= NODE_ON_BEST_BRANCH;}
    void setNodeOnNextBestBranch()               {m_statusBits |= NODE_ON_NEXT_BEST_BRANCH;}
    void setBiLyrs2MainBrch(int biLyrs);

    // Methods to return information
    // Return pointer to the parent node
    const TkrVecNode*                   getParentNode()      const {return m_parent;}
    // Return pointer to the associated link
    const TkrVecPointsLink*             getAssociatedLink()  const {return m_associatedLink;}
    // Return the rms angle
    const double                        getRmsAngle()        const;
    // Return the depth of tree to this node
    const int                           getDepth()           const {return m_depth;}
    // Number of leaves in the tree
    const int                           getNumLeaves()       const {return m_leaves;}
    // number of branches in the tree
    const int                           getNumBranches()     const {return m_branches;}
    // Number of nodes in the tree
    const int                           getNumNodesInTree()  const;
    // Return the memory used by the tree
    const int                           getMemUsedInTree()   const;
    // Get the tree starting layer
    const unsigned char                 getTreeStartLayer()  const;
    // Get the current bilayer
    const unsigned char                 getCurrentBiLayer()  const;
    // Get the current bilayer
    const unsigned char                 getBiLyrs2MainBrch() const;
    // Get the Tree ID
    const unsigned char                 getTreeId()          const;
    // Get the number of bilayers from start to this node
    const unsigned char                 getNumBiLayers()     const {return getTreeStartLayer() - getCurrentBiLayer() + 1;}
    // Check status bits
    const bool                          isNodeSacred()       const {return (m_statusBits & NODE_IS_SACRED)           != 0;}
    const bool                          canNodeBeShared()    const {return (m_statusBits & NODE_CAN_BE_SHARED)       != 0;}
    const bool                          isOnBestBranch()     const {return (m_statusBits & NODE_ON_BEST_BRANCH)      != 0;}
    const bool                          isOnNextBestBranch() const {return (m_statusBits & NODE_ON_NEXT_BEST_BRANCH) != 0;}
    // Get the number of BiLayers from this node along best branch
    const int                           getBestNumBiLayers() const {return m_bestNumBiLayers;}
    // Get the rms deviations from this node along best branch
    const double                        getBestRmsAngle()    const; // {return m_bestRmsAngle;}

    // Define operator to facilitate sorting
    const bool operator<(const TkrVecNode& right) const
    {
        // kludge in a test here... this needs work...
        // The tree that starts higher gets precedence
        if      (this->getTreeStartLayer() > right.getTreeStartLayer()) return true;
        else if (this->getTreeStartLayer() < right.getTreeStartLayer()) return false;
        // Both angles non-zero
        else if (m_rmsAngleSum > 0. && right.getRmsAngleSum() > 0.)
        {
            double myPenalty     = 1.; //m_associatedLink->skipsLayers() ? m_associatedLink->skip1Layer() ? 1.25 : 1.50 : 1.0;
            double myRmsAngle    = myPenalty * m_rmsAngleSum / m_numAnglesInSum;
            double rightPenalty  = 1.; //right.getAssociatedLink()->skipsLayers() ? right.getAssociatedLink()->skip1Layer() ? 1.25 : 1.50 : 1.0;
            double rightRmsAngle = rightPenalty * right.getRmsAngle();

            return myRmsAngle < right.getRmsAngle();
        }
        else if (right.getRmsAngleSum() > 0.) return false;

        return true;
    }

    // Direct access to angle sums 
    const double getRmsAngleSum()    const {return m_rmsAngleSum;}
    const int    getNumAnglesInSum() const {return m_numAnglesInSum;}

    // "Proper" removal of a daughter node
    bool removeDaughter(Event::TkrVecNode* node, bool keepBest = true);

private:
    // For resetting the "best" parameters by going up the parent tree
    void resetBestParams();

    // Pointer to the link we are associated to
    const TkrVecPointsLink*      m_associatedLink;
    // Pointer to our parent node
    TkrVecNode*                  m_parent;
    // Bit mask to contain status items
    // Currently: bits  0- 4 contain the tree starting bilayer
    //            bits  5- 9 contain the current bilayer
    //            bits 10-14 contian the number of bilayers to main branch
    //            bits 16-23 contain the tree ID
    //            bits 24-31 represent status bits:
    //            0x80000000 is the "this node is sacred" bit
    unsigned int                 m_statusBits;
    // Following give information TO (and including) this node
    double                       m_rmsAngleSum;
    int                          m_numAnglesInSum;
    // Following give information FROM (and including) this node
    int                          m_leaves;             // Number of leaves below this node
    int                          m_branches;           // Number of branches (total) below this node
    int                          m_depth;              // "Depth" of tree below this node
    int                          m_bestNumBiLayers;    // Number of bilayers below this node for best branch
    double                       m_bestRmsAngle;       // RMS deviation of links below for best branch
};


inline TkrVecNode::TkrVecNode(TkrVecNode* parent, const TkrVecPointsLink* associatedLink) : 
                              m_associatedLink(associatedLink),
                              m_parent(parent),
                              m_statusBits(0),
                              m_rmsAngleSum(0.),
                              m_numAnglesInSum(1),
                              m_leaves(0),
                              m_branches(0),
                              m_depth(1),
                              m_bestNumBiLayers(0),
                              m_bestRmsAngle(0.)
{
    // If we have a parent node then add ourselves as daughter
    // and update our angle info as the base
    if (parent)
    {
        // Under new scheme we can "self insert" outselves into the parent's daughter list since we
        // don't have a corresponding TkrVecNodeInfo object...
//        m_parent->addDaughter(this);

        m_rmsAngleSum    = parent->getRmsAngleSum();
        m_numAnglesInSum = parent->getNumAnglesInSum();

        if (parent->getTreeStartLayer())
        {
            m_statusBits = (m_statusBits & ~START_BILAYER_BITS) 
                         | (parent->getTreeStartLayer() & START_BILAYER_BITS);
        }
        else if (associatedLink)
        {
            int biLayer = associatedLink->getFirstVecPoint()->getLayer();

            m_statusBits = (m_statusBits & START_BILAYER_BITS) | (biLayer & START_BILAYER_BITS);
        }

        setTreeId(parent->getTreeId());
    }

    // If we have an associated link then update the angle information
    if (associatedLink && parent->getAssociatedLink())
    {
        TkrVecPointsLink* parentLink = const_cast<Event::TkrVecPointsLink*>(parent->getAssociatedLink());

        double linkAngle = parentLink->angleToNextLink(*associatedLink);
//        double linkAngle = parentLink->distToNextLink(*associatedLink);

        m_rmsAngleSum += linkAngle * linkAngle;
        m_numAnglesInSum++;

        if (associatedLink->skipsLayers())
        {
            m_rmsAngleSum += linkAngle * linkAngle;

            if (associatedLink->skip2Layer())
            {
                m_rmsAngleSum += linkAngle * linkAngle;
            }
        }
    }

    // Set the current bilayer
    if (associatedLink)
    {
        int biLayer = associatedLink->getFirstVecPoint()->getLayer();

        m_statusBits = (m_statusBits & ~CURRENT_BILAYER_BITS) | ((biLayer << 5) & CURRENT_BILAYER_BITS);
    }
}

inline TkrVecNode::~TkrVecNode()
{
    // Set the contained object parent status to zero
    setParent(0);

    // We go through and delete our daughters. Note that the daughters will "remove" themselves
    // from this list in their destructor, hence the odd looking looping mechanism to try to
    // insure we have a valid iterator always
    TkrVecNodeSet::iterator nodeItr = begin(); 
    while(nodeItr != end())
    {
        // Dereference the node pointer
        TkrVecNode* node = *nodeItr++;

        // Delete it
        delete node;
    }

    // Its either now or in the list destructor...
    clear();

    // If we have a parent, then remove ourselves as a daughter of mom or pop
    // Well, we can self-remove under new system either. hmmm
    if (m_parent) m_parent->removeDaughter(this);

    return;
}

inline void TkrVecNode::setRmsAngleInfo(double rmsAngSum, int num)   
{
    m_rmsAngleSum    = rmsAngSum;
    m_numAnglesInSum = num;

    return;
}

inline const double TkrVecNode::getRmsAngle() const
{
    double rmsAngle = sqrt(m_rmsAngleSum / m_numAnglesInSum);

    return rmsAngle;
}
    
inline const double TkrVecNode::getBestRmsAngle() const
{
    double rmsAngle = m_bestRmsAngle; // * double(getNumBiLayers()) / double(m_depth);

    return rmsAngle;
}

//inline const int TkrVecNode::getDepth() const
//{
//    int depthBelowMe = 0;
//
//    if (!empty())
//    {
//        // Insertion of daughters is sorted, first node should be the "deepest"
//        depthBelowMe = (*begin())->getDepth();
//    }
//
//    return depthBelowMe + 1;
//}

//inline const int TkrVecNode::getNumLeaves() const
//{
//    // If we have no daughters then we are the leaf and need to be counted
//    int leafCount = 1;
//
//    if (!empty())
//    {
//        leafCount = 0;
//
//        for(TkrVecNodeSet::const_iterator nodeItr = begin(); nodeItr != end(); nodeItr++)
//        {
//            leafCount += (*nodeItr)->getNumLeaves();
//        }
//    }
//
//    return leafCount;
//}

inline const int TkrVecNode::getNumNodesInTree() const
{
    // We always count ourself
    int nodeCount = m_associatedLink ? 1 : 0;

    if (!empty())
    {
        for(TkrVecNodeSet::const_iterator nodeItr = begin(); nodeItr != end(); nodeItr++)
        {
            nodeCount += (*nodeItr)->getNumNodesInTree();
        }
    }

    return nodeCount;
}

inline const int TkrVecNode::getMemUsedInTree() const
{
    // We always count ourself
    int memUsed = sizeof(*this) + size() * sizeof(Event::TkrVecNode*);

    if (!empty())
    {
        for(TkrVecNodeSet::const_iterator nodeItr = begin(); nodeItr != end(); nodeItr++)
        {
            memUsed += (*nodeItr)->getMemUsedInTree();
        }
    }

    return memUsed;
}
    
inline const unsigned char TkrVecNode::getTreeStartLayer() const
{
    unsigned char startLayer = m_statusBits & START_BILAYER_BITS;

    return startLayer;
}

inline const unsigned char TkrVecNode::getCurrentBiLayer()  const
{
    unsigned char currentLayer = (m_statusBits & CURRENT_BILAYER_BITS) >> 5;

    return currentLayer;
}

inline const unsigned char TkrVecNode::getBiLyrs2MainBrch()  const
{
    unsigned char currentLayer = (m_statusBits & TO_MAIN_BRANCH_BITS) >> 10;

    return currentLayer;
}

inline const unsigned char TkrVecNode::getTreeId() const
{
    unsigned char treeId = (m_statusBits & TREE_ID_BITS) >> 16;

    return treeId;
}

inline void TkrVecNode::setTreeId(int treeId)
{
    m_statusBits = (m_statusBits & ~TREE_ID_BITS) | ((treeId << 16) & TREE_ID_BITS);
}

inline void TkrVecNode::setBiLyrs2MainBrch(int biLayersToMainBranch)
{
    m_statusBits = (m_statusBits & ~TO_MAIN_BRANCH_BITS) | ((biLayersToMainBranch << 10) & TO_MAIN_BRANCH_BITS);
}

inline bool TkrVecNode::removeDaughter(Event::TkrVecNode* node, bool keepBest)
{
    bool removed = false;

//    Event::TkrVecNodeSet::iterator nodeIter = find(node);
    Event::TkrVecNodeSet::iterator nodeIter = std::find(begin(), end(), node);

    // Protection against attempting to remove a node which is not a daughter
    if (nodeIter != end())
    {
        // Special case is that node to be removed is the first node (unlikely?)
        bool resetBest = keepBest && nodeIter == begin() ? true : false;

        erase(nodeIter);

        // If that was the "best" node then reset the best params
        if (resetBest) resetBestParams();
        
        // If no daughters then we are the new end of the line
        if (empty())
        {
            m_leaves          = 1;
            m_branches        = 0;
        }
        // Otherwise we need to update our information based on remaining daughters
        else
        {
            // Reset the number of leaves, etc.
            m_leaves   -= node->getNumLeaves();
            m_branches -= node->getNumBranches();
        }

        removed = true;
    }

    return removed;
}

inline void TkrVecNode::resetBestParams()
{
    // If our daughter list is not empty then reset the "best" params to match the first
    // daughter in the list. 
    // **** NOTE: I need to be fixed to prevent case where caller is no longer first in list
    if (!empty())
    {
        TkrVecNode* node = *begin();

        m_depth           = node->getDepth() + 1;
        m_bestNumBiLayers = node->getBestNumBiLayers();
        m_bestRmsAngle    = node->getBestRmsAngle();
    }
    else
    {
        m_depth           = 1;
        m_bestNumBiLayers = getTreeStartLayer() - getCurrentBiLayer() + 1;
        m_bestRmsAngle    = getRmsAngle();
    }

    // Reset on up the chain of command...
    if (m_parent) m_parent->resetBestParams();

    return;
}


inline const bool TkrVecNodesComparator::operator()(const TkrVecNode* left, const TkrVecNode* right) const
{
    // Most number of bilayers wins (longest)
    if      (left->getBestNumBiLayers() > right->getBestNumBiLayers()) return true;
    else if (left->getBestNumBiLayers() < right->getBestNumBiLayers()) return false;

    // Check special case of stubs starting with skipping layer links
    if (left->getNumAnglesInSum() == 2 || right->getNumAnglesInSum() == 2)
    {
        if      (left->getDepth() < right->getDepth()) return false;
        else if (left->getDepth() > right->getDepth()) return true;
    }

    // Last check is to take the branch that is "straightest". 
    // Use the scaled rms angle to determine straightest...
    double leftRmsAngle  = left->getBestRmsAngle() * double(left->getNumBiLayers()) / double(left->getDepth());
    double rightRmsAngle = right->getBestRmsAngle() * double(right->getNumBiLayers()) / double(right->getDepth());
    
    //if (left->getBestRmsAngle() < right->getBestRmsAngle()) return true;
    if (leftRmsAngle < rightRmsAngle) return true;

    return false;
}


// These typedefs useful for point/link association
typedef std::list<TkrVecNode>                      TkrVecNodeList;
typedef std::list<TkrVecNode*>                     TkrVecNodePtrList;

// Typedefs for gaudi container for these objects
typedef ObjectVector<TkrVecNode>                   TkrVecNodeCol;
typedef TkrVecNodeCol::iterator                    TkrVecNodeColPtr;
typedef TkrVecNodeCol::const_iterator              TkrVecNodeColConPtr;

typedef RelTable<TkrVecPoint, TkrVecNode>          TkrVecPointToNodesTab;
typedef Relation<TkrVecPoint, TkrVecNode>          TkrVecPointToNodesRel;
typedef RelationList<TkrVecPoint, TkrVecNode>      TkrVecPointToNodesTabList;

typedef RelTable<const TkrCluster, TkrVecNode>     TkrClusterToNodesTab;
typedef Relation<const TkrCluster, TkrVecNode>     TkrClusterToNodesRel;
typedef RelationList<const TkrCluster, TkrVecNode> TkrClusterToNodesTabList;

}; // Namespace

#endif

