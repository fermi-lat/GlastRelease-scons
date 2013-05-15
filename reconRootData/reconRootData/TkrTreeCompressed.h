/** @file TkrTreeCompressed.h

* @class TkrTreeCompressed
*
* @brief This is the root version of the TkrTreeCompressed object. The TkrTreeCompressed structure is the top
*        level object returned by the Tree Based Tracking pattern recognition
*
* last modified 02/01/2013
*
* @authors Usher
*
* $Header: 
*/

#ifndef ROOT__TkrTreeCompressed_H
#define ROOT__TkrTreeCompressed_H

#include "TObjArray.h"
#include "reconRootData/TkrVecNodesCompressed.h"
#include "reconRootData/TkrFilterParams.h"
#include "reconRootData/TkrTrack.h"
#include <vector>
#ifndef R__GLOBALSTL
#ifndef WIN32
using std::vector;
#else
using namespace std;
#endif
#endif

class TkrTreeCompressed: public TObjArray
{
public:
    // constructors - null case
    TkrTreeCompressed();

    // constructurs - have the information case
//    TkrTreeCompressed(TkrTreeCompressed* parent, const TkrVecPointsLink* associatedLink);

    // destructor
    virtual ~TkrTreeCompressed() {Clear();}

    // Root mandated methods
    void Clear(Option_t *option="");
    void Fake( Int_t ievent, UInt_t rank, Float_t randNum ) ; // for tests
    Bool_t CompareInRange( const TkrTreeCompressed &, const std::string & name = "" ) const ; // for tests
    void Print(Option_t *option="") const;

    // Initializer
    void initializeInfo(TkrVecNodeCompressed* head, 
                        TkrVecNodeCompressed* best, 
                        TkrVecNodeCompressed* second, 
                        TkrFilterParams*      filter,
                        Double_t              bestBranchAngToAxis,
                        Double_t              axisSeededAngToAxis);

    // Methods to return information
    // Return pointer to the head node for this tree
    const TkrVecNodeCompressed* getHeadNode()              const {return m_headNode;}
    // Return pointer to the best leaf, corresponding to the best track
    const TkrVecNodeCompressed* getBestLeaf()              const {return m_bestLeaf;}
    // Return pointer to the second leaf, corresponding to second track (if there)
    const TkrVecNodeCompressed* getSecondLeaf()            const {return m_secondLeaf;}
    // Return pointer to the sibling map
//    const TkrNodeSiblingMap* getSiblingMap()            const {return m_siblingMap;}
    // Return pointer to the axis parameters
    const TkrFilterParams*      getAxisParams()            const {return m_axisParams;}
    // Return pointer to the track for this tree
    const TkrTrack*             getBestTrack()             const {return GetEntries() > 0 ? dynamic_cast<TkrTrack*>(At(0)) : 0; }
    // Return angle best branch track makes with tree axis
    const Double_t              getBestBranchAngleToAxis() const {return m_bestBranchAngleToAxis;}
    // Return angle axis seeded track makes with tree axis
    const Double_t              getAxisSeededAngleToAxis() const {return m_axisSeededAngleToAxis;}

    /// Override the "new" operator in order to use the ReconObjectManager pool management
    void* operator new(size_t size);

    void* operator new(size_t size, void* vp);

    /// Override the "delete" operator in order to use the ReconObjectManager pool management
    void  operator delete(void* p);

private:
    // Pointer to the head node in the Tree
    TkrVecNodeCompressed* m_headNode;

    // Pointer to the "best" leaf corresponding to "best" track (if not composite)
    TkrVecNodeCompressed* m_bestLeaf;

    // Pointer to the second leaf, corresponding to the second track (if exists)
    TkrVecNodeCompressed* m_secondLeaf;

    // Pointer to the Node sibling map
//    TkrNodeSiblingMap* m_siblingMap;

    // Pointer to the tree axis parameters 
    TkrFilterParams*      m_axisParams;

    // Angle best leaf track to axis
    Double32_t            m_bestBranchAngleToAxis;

    // Angle axis seeded track to axis
    Double32_t            m_axisSeededAngleToAxis;


    ClassDef(TkrTreeCompressed,1)
};


#endif // __TkrTreeCompressed_H
