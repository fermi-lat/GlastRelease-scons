/** @file TkrVecNodes.h

* @class TkrVecNode
*
* @brief This is the root version of the TkrVecNode object. The purpose of the TkrVecNode object is
*        to represent a 3D "point" in the Tracker definied by the crossing of two TkrClusters
*
* last modified 01/29/2013
*
* @authors Usher
*
* $Header: 
*/

#ifndef ROOT__TkrVecNodes_H
#define ROOT__TkrVecNodes_H

#include "TObjArray.h"
#include "reconRootData/TkrVecPointsLink.h"
#include <vector>
#ifndef R__GLOBALSTL
#ifndef WIN32
using std::vector;
#else
using namespace std;
#endif
#endif

class TkrVecNode: public TObjArray
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

    // constructors - null case
    TkrVecNode() {Clear();}

    // destructor
    virtual ~TkrVecNode() {Clear();}

    // Root mandated methods
    void Clear(Option_t *option="");
    void Fake( Int_t ievent, UInt_t rank, Float_t randNum ) ; // for tests
    Bool_t CompareInRange( const TkrVecNode &, const std::string & name = "" ) const ; // for tests
    void Print(Option_t *option="") const;

    // Initializer
    void initializeInfo(TkrVecNode*             parent,
                        const TkrVecPointsLink* associatedLink, 
                        UInt_t                  statusBits,
                        Double_t                rmsAngleSum,
                        Int_t                   numAnglesInSum,
                        Int_t                   leaves,
                        Int_t                   branches,
                        Int_t                   depth,
                        Int_t                   bestNumBiLayers,
                        Double_t                bestRmsAngle);

    // Methods to manipulate the data members
    // Method to update the rms angle going TO this node
    void setRmsAngleInfo(double rms, int num);
    // Set the number of leaves below this node
    void setNumLeaves(Int_t numLeaves)           {m_leaves          = numLeaves;}
    // Set the number of branches below this nodde
    void setNumBranches(Int_t numBranches)       {m_branches        = numBranches;}
    // Set the depth below this node
    void setDepth(Int_t depth)                   {m_depth           = depth;}
    // Set the number of bilayers below this node for best branch
    void setBestNumBiLayers(Int_t numBiLayers)   {m_bestNumBiLayers = numBiLayers;}
    // Set the RMS deviation below this node for best branch
    void setBestRmsAngle(Double_t rmsAngle)      {m_bestRmsAngle    = rmsAngle;}
    // Set the Tree ID #
    void setTreeId(Int_t treeId);
    // Set status bits
    void setNodeSacred()                         {m_statusBits |= NODE_IS_SACRED;}
    void setNodeCanBeShared()                    {m_statusBits |= NODE_CAN_BE_SHARED;}
    void setNodeOnBestBranch()                   {m_statusBits |= NODE_ON_BEST_BRANCH;}
    void setNodeOnNextBestBranch()               {m_statusBits |= NODE_ON_NEXT_BEST_BRANCH;}
    void setTreeStartLayer(Int_t startLayer);
    void setBiLyrs2MainBrch(Int_t biLyrs);

    // Methods to return information
    // Return pointer to the parent node
    const TkrVecNode*                   getParentNode()           const {return m_parent;}
    // Return pointer to the associated link
    const TkrVecPointsLink*             getAssociatedLink()       const {return m_associatedLink;}
    // Return the rms angle
    const Double_t                      getRmsAngle()             const;
    // Return the depth of tree to this node
    const Int_t                         getDepth()                const {return m_depth;}
    // Number of leaves in the tree
    const Int_t                         getNumLeaves()            const {return m_leaves;}
    // Get the number of branches
    const Int_t                         getNumBranches()          const {return m_branches;}
    // Get the tree starting layer
    const unsigned char                 getTreeStartLayer()       const;
    // Get the current bilayer
    const unsigned char                 getCurrentBiLayer()       const;
    // Get the current bilayer
    const unsigned char                 getBiLyrs2MainBrch()      const;
    // Get the Tree ID
    const unsigned char                 getTreeId()               const;
    // Get the number of bilayers from start to this node
    const unsigned char                 getNumBiLayers()          const {return getTreeStartLayer() - getCurrentBiLayer() + 2;}
    // Check status bits
    const bool                          isNodeSacred()            const {return (m_statusBits & NODE_IS_SACRED)           != 0;}
    const bool                          canNodeBeShared()         const {return (m_statusBits & NODE_CAN_BE_SHARED)       != 0;}
    const bool                          isOnBestBranch()          const {return (m_statusBits & NODE_ON_BEST_BRANCH)      != 0;}
    const bool                          isOnNextBestBranch()      const {return (m_statusBits & NODE_ON_NEXT_BEST_BRANCH) != 0;}
    // Return the full Status Bits
    const UInt_t                        getStatusBits()           const {return m_statusBits;}
    // Get the number of BiLayers from this node along best branch
    const Int_t                         getBestNumBiLayers()      const {return m_bestNumBiLayers;}
    // Get the rms deviations from this node along best branch
    const Double_t                      getBestRmsAngle()         const; // {return m_bestRmsAngle;}

    // Direct access to angle sums 
    const Double_t                      getRmsAngleSum()          const {return m_rmsAngleSum;}
    const Int_t                         getNumAnglesInSum()       const {return m_numAnglesInSum;}

    /// Override the "new" operator in order to use the ReconObjectManager pool management
    void* operator new(size_t size);

    void* operator new(size_t size, void* vp);

    /// Override the "delete" operator in order to use the ReconObjectManager pool management
    void  operator delete(void* p);

private:
    // Pointer to the link we are associated to
    const TkrVecPointsLink* m_associatedLink;
    // Pointer to our parent node
    TkrVecNode*             m_parent;
    // Bit mask to contain status items
    // Currently: bits  0- 4 contain the tree starting bilayer
    //            bits  5- 9 contain the current bilayer
    //            bits 10-14 contian the number of bilayers to main branch
    //            bits 16-23 contain the tree ID
    //            bits 24-31 represent status bits:
    //            0x80000000 is the "this node is sacred" bit
    UInt_t                  m_statusBits;
    // Following give information TO (and including) this node
    Double32_t              m_rmsAngleSum;
    Int_t                   m_numAnglesInSum;
    // Following give information FROM (and including) this node
    Int_t                   m_leaves;             // Number of leaves below this node
    Int_t                   m_branches;           // Number of branches (total) below this node
    Int_t                   m_depth;              // "Depth" of tree below this node
    Int_t                   m_bestNumBiLayers;    // Number of bilayers below this node for best branch
    Double32_t              m_bestRmsAngle;       // RMS deviation of links below for best branch

    ClassDef(TkrVecNode,1)
};


#endif // __TkrVecNode_H
