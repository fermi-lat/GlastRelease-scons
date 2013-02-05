#ifndef ROOT_TKRRECON_H
#define ROOT_TKRRECON_H

#include "TObject.h"
#include "TObjArray.h"

#include "TkrCluster.h"
#include "TkrTree.h"
#include "TkrTrack.h"
#include "TkrVertex.h"
#include "TkrTruncationData.h"
#include "TkrDiagnostics.h"
#include "TkrFilterParams.h"
#include "TkrEventParams.h"
#include "TreeClusterRelation.h"
/** 
 * @class TkrRecon
 * @brief Root container for the reconstructed Tracker information 
 * 
 * This class contains:
 * -# Collection of Clusters
 * -# Collection of 
 * -# Collectoin of Reconstructed Tracks
 * -# Collection Vertices
 * 
 * @author Tracy Usher
 *
 * $Header$
 */

class TkrRecon : public TObject
{

public:
    TkrRecon();

    virtual ~TkrRecon();

    void initialize();
    
    void Clear(Option_t *option="");
    void Fake( Int_t ievent, Float_t randNum ) ; // for tests
    Bool_t CompareInRange( TkrRecon &, const std::string & name = "" ) ; // for tests
    void Print(Option_t *option="") const;

    /// provide access to the list of clusters    
    TObjArray *getClusterCol() { return m_clusterCol; };
    void addCluster(TkrCluster *cluster) { m_clusterCol->Add(cluster); };

    TObjArray *getTrackCol() { return m_trackCol; };
    void addTrack(TObject* track) { m_trackCol->Add(track); };
    /// ADW: Added cosmic-ray track collection (eventually make a map of collections?)
    TObjArray *getCRTrackCol() { return m_crTrackCol; };
    void addCRTrack(TObject* track) { m_crTrackCol->Add(track); };

    TObjArray *getVertexCol() { return m_vertexCol; };
    void addVertex(TkrVertex* vertex) { m_vertexCol->Add(vertex); };

    TObjArray *getTruncationDataCol() { return m_truncationDataCol; };
    void addTruncationData(TkrTruncationData* trdata) { m_truncationDataCol->Add(trdata); };

    TObjArray* getTkrVecPointCol()    {return m_tkrVecPointCol;}
    void addTkrVecPoint(TkrVecPoint* vecPoint) {m_tkrVecPointCol->Add(vecPoint);}

    TObjArray* getTkrVecPointsLinkCol()    {return m_tkrVecPointsLinkCol;}
    void addTkrVecPointsLink(TkrVecPointsLink* vecPoint) {m_tkrVecPointsLinkCol->Add(vecPoint);}

    TObjArray* getTkrVecNodesCol()    {return m_tkrVecNodesCol;}
    void addTkrVecNode(TkrVecNode* vecNode) {m_tkrVecNodesCol->Add(vecNode);}

    TObjArray* getTkrFilterParamsCol()    {return m_tkrFilterParamsCol;}
    void addTkrFilterParams(TkrFilterParams* filterParams) {m_tkrFilterParamsCol->Add(filterParams);}

    TkrEventParams* getTkrEventParams()   {return m_tkrEventParams;}
    void addTkrEventParams(TkrEventParams* eventParams) {m_tkrEventParams = eventParams;}

    TObjArray* getTkrTreeCol()    {return m_tkrTreeCol;}
    void addTkrTree(TkrTree* tree) {m_tkrTreeCol->Add(tree);}

    TObjArray* getTreeClusterRelationCol()    {return m_treeClusterRelationCol;}
    void addTreeClusterRelation(TreeClusterRelation* treeClusterRelation) {m_treeClusterRelationCol->Add(treeClusterRelation);}

    // Added for macros & interactive access - EAC
    Int_t nCluster() const { return m_clusterCol->GetEntries() ; }    
    Int_t nTrack() const { return m_trackCol->GetEntries() ; }
    Int_t nCRTrack() const { return m_crTrackCol->GetEntries() ; }
    Int_t nVertex() const { return m_vertexCol->GetEntries() ; }
    Int_t nTruncationData() const { return m_truncationDataCol->GetEntries() ; }
    Int_t nTkrVecPoints() const { return m_tkrVecPointCol->GetEntries();}
    Int_t nTkrVecPointsLinks() const { return m_tkrVecPointsLinkCol->GetEntries();}
    Int_t nTkrVecNodes() const {return m_tkrVecNodesCol->GetEntries();}
    Int_t nTkrTrees() const {return m_tkrTreeCol->GetEntries();}
    Int_t nTkrFilterParams() const {return m_tkrFilterParamsCol->GetEntries();}
    Int_t nTkrEventParams() const {return m_tkrEventParams ? 1 : 0;}
    Int_t nTreeClusterRelations() const {return m_treeClusterRelationCol->GetEntries();}

    TkrCluster* getCluster(Int_t i) { return i < nCluster() ? dynamic_cast<TkrCluster*>(m_clusterCol->At(i)) : 0; };
    TkrTrack* getTrack(Int_t i) { return i < nTrack() ? dynamic_cast<TkrTrack*>(m_trackCol->At(i)) : 0; };
    TkrTrack* getCRTrack(Int_t i) { return i < nCRTrack() ? dynamic_cast<TkrTrack*>(m_crTrackCol->At(i)) : 0; };
    TkrVertex* getVertex(Int_t i) { return i < nVertex() ? dynamic_cast<TkrVertex*>(m_vertexCol->At(i)) : 0; };
    TkrTruncationData* getTruncationData(Int_t i) { return i < nTruncationData() ? dynamic_cast<TkrTruncationData*>(m_truncationDataCol->At(i)) : 0; };
    TkrVecPoint* getTkrVecPoint(Int_t i) { return i < nTkrVecPoints() ? dynamic_cast<TkrVecPoint*>(m_tkrVecPointCol->At(i)) : 0;}
    TkrVecPointsLink* getTkrVecPointsLink(Int_t i) { return i < nTkrVecPointsLinks() ? dynamic_cast<TkrVecPointsLink*>(m_tkrVecPointsLinkCol->At(i)) : 0;}
    TkrVecNode* getTkrVecNode(Int_t i) { return i < nTkrVecNodes() ? dynamic_cast<TkrVecNode*>(m_tkrVecNodesCol->At(i)) : 0;}
    TkrTree* getTkrTree(Int_t i) { return i < nTkrTrees() ? dynamic_cast<TkrTree*>(m_tkrTreeCol->At(i)) : 0;}
    TkrFilterParams* getTkrFilterParams(Int_t i) { return i < nTkrFilterParams() ? dynamic_cast<TkrFilterParams*>(m_tkrFilterParamsCol->At(i)) : 0;}
    TkrEventParams* getTkrEventParams(Int_t i) { return m_tkrEventParams;}
    TreeClusterRelation* getTreeClusterRelation(Int_t i) {return i < nTreeClusterRelations() ? dynamic_cast<TreeClusterRelation*>(m_treeClusterRelationCol->At(i)) : 0;}
    

    const TkrDiagnostics* getDiagnostics() { return m_diagnostics;};
    void addDiagnostics(TkrDiagnostics* diagnostics) { m_diagnostics = diagnostics;};

private:
    /// collection of TkrSiCluster
    TObjArray *m_clusterCol;
    /// collection of TkrTrack
    TObjArray *m_trackCol;
    /// collection of TkrVertex
    TObjArray *m_vertexCol;
     /// Diagnostics
    TkrDiagnostics*   m_diagnostics;
    /// collection of TkrTruncationData
    TObjArray *m_truncationDataCol;
    /// ADW: cosmic-ray collection of TkrTrack
    TObjArray *m_crTrackCol;
    /// collection of TkrVecPoints
    TObjArray*      m_tkrVecPointCol;
    /// collection of TkrVecPointsLink
    TObjArray*      m_tkrVecPointsLinkCol;
    /// collection of TkrVecNodes
    TObjArray*      m_tkrVecNodesCol;
    /// collection of TkrFilterParams
    TObjArray*      m_tkrFilterParamsCol;
    /// collection of TkrTrees
    TObjArray*      m_tkrTreeCol;
    /// Pointer to the TkrEventParams object (only one per event)
    TkrEventParams* m_tkrEventParams;
    /// Pointer to the TreeClusterRelations
    TObjArray*      m_treeClusterRelationCol;


    static TkrCluster *keepCluster[16384]; //!
    static TkrTrack *keepTrack[10000]; //!
    static TkrVertex *keepVertex[10000]; //!
    // HMK-v12r7 static TkrTruncationData *keepTruncationData[10000]; //!
    static Int_t indCluster, indTrack, indVertex; //!        
    //, indTruncationData; 

    ClassDef(TkrRecon,10)
};

#endif
