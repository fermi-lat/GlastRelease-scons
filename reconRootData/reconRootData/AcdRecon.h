#ifndef ROOT_ACDRECON_H
#define ROOT_ACDRECON_H

#include "TObject.h"
#include "commonRootData/idents/AcdId.h"

#include "AcdTkrIntersection.h"
#include "TClonesArray.h"

#include <vector>
#include "TClonesArray.h"

#include "AcdTkrIntersection.h"

#ifndef R__GLOBALSTL
#ifndef WIN32
using std::vector;
#else
using namespace std;
#endif
#endif

/** 
* @class AcdRecon
* @brief Root container for the reconstructed ACD information 
* 
* This class contains:
* -# Minimum DOCA (Distance of Closest Approach) computed using all tracks and
*    ACD detectors.  DOCA is computed using the center of the ACD tiles.
* -# AcdId of the tile associated with the minimum DOCA.
* -# Minimum Active Distance computed using all tracks and ACD detectors.
*    Computed using the edge of ACD tiles.
* -# AcdId of the tile associated with the maximum Active Distance
* -# Tile Count the number of ACD tiles which were above veto threshold.
* -# Collection of minimum DOCA values organized by geometry:  top, side rows
* -# Collection of minimum Active Distance values organized by geometry:  top
*    and side rows.
* -# Total reconstructed energy detected in ACD detector elements.
* -# MC Energy in MeV deposited in the ribbons
* -# Ribbon Count
*
* @author Heather Kelly
*
* $Header$
*/

class AcdRecon : public TObject 
{
public:
    AcdRecon();
 
    AcdRecon(Double_t e, Double_t ribbonE, Int_t count, Int_t ribbonCount, 
            Double_t gDoca, Double_t doca, 
            const AcdId &minDocaId, Double_t actDist, const AcdId &maxActDistId,
            const std::vector<Double_t> &rowDoca,
            const std::vector<Double_t> &rowActDist,
            const std::vector<AcdId> &idCol, 
            const std::vector<Double_t> &energyCol)       
            : m_totEnergy(e),
            m_gammaDoca(gDoca),
            m_doca(doca),
            m_minDocaId(minDocaId),
            m_actDist(actDist),
            m_maxActDistId(maxActDistId),
            m_tileCount(count),
            m_rowDocaCol(rowDoca),
            m_rowActDistCol(rowActDist),
            m_idCol(idCol),
            m_energyCol(energyCol),
            m_totRibbonEnergy(ribbonE),
            m_ribbonCount(ribbonCount)
        {};
     
    virtual ~AcdRecon();
    
    void initialize(Double_t e, Double_t ribbonE, Int_t count, 
        Int_t ribbonCount, Double_t gDoca, Double_t doca, 
        const AcdId &minDocaId, Double_t actDist, 
        const AcdId &maxActDistId, Double32_t ribbonActDist,
        const AcdId &ribbonActDistId, const std::vector<Double_t> &rowDoca, 
        const std::vector<Double_t> &rowActDist, 
        const std::vector<AcdId> &idCol,
        const std::vector<Double_t> &energyCol);

    /// overload initialize for interactive root where we cannot use std::vector
    void initialize(Double_t e, Double_t ribbonE, Int_t count, 
                    Int_t ribbonCount, Double_t gDoca, Double_t doca, 
                    const AcdId &minDocaId,
                    Double_t actDist, const AcdId &maxActDistId);
    
    void Clear(Option_t *option="");
    
    void Print(Option_t *option="") const;
    
    inline const Double_t getEnergy() const { return m_totEnergy; };
    inline const Double_t getRibbonEnergy() const { return m_totRibbonEnergy; };
    inline const Int_t getTileCount() const { return m_tileCount; };
    inline const Int_t getRibbonCount() const { return m_ribbonCount; };
    inline const Double_t getGammaDoca() const { return m_gammaDoca; };
    inline const Double_t getDoca() const { return m_doca; };
    inline const Double_t getActiveDist() const { return m_actDist; };
    inline const Double32_t getRibbonActiveDist() const { return m_ribbonActDist; };
    inline const AcdId& getMinDocaId() const { return m_minDocaId; };
    inline const AcdId& getMaxActDistId() const { return m_maxActDistId; };
    inline const AcdId& getRibbonActDistId() const { return m_ribbonActDistId; }

    inline const std::vector<Double_t>& getRowDocaCol() const { return m_rowDocaCol; };
    inline const Double_t getRowDoca(UInt_t i) const { 
        return ((i < m_rowDocaCol.size()) ? m_rowDocaCol[i] : -1.); };
    inline void addRowDoca(Double_t val) { m_rowDocaCol.push_back(val); };
    inline const std::vector<Double_t>& getRowActDistCol() const { return m_rowActDistCol; };
    inline const Double_t getRowActDist(UInt_t i) const {
        return ( (i < m_rowActDistCol.size()) ? m_rowActDistCol[i] : -1.); };
    inline void addRowActDist(Double_t val) { m_rowActDistCol.push_back(val); };
    inline const std::vector<Double_t>& getEnergyCol() const { return m_energyCol; };
    inline const Double_t getEnergy(UInt_t i) const { 
           return ( (i < m_energyCol.size()) ? m_energyCol[i] : -1.); };
    inline void addEnergy(Double_t e) { m_energyCol.push_back(e); };
    inline const std::vector<AcdId>& getIdCol() const { return m_idCol; };
    inline const AcdId* getId(UInt_t i) const { 
           return ( (i < m_idCol.size()) ? &m_idCol[i] : 0); };
    inline void addId(const AcdId &id) { m_idCol.push_back(id); };
    Double_t getEnergy(const AcdId &id) const;

    // for ROI stuff
    Bool_t hasHit(UInt_t acdId) const;
    inline UInt_t nId() { return m_idCol.size(); };
  
    // does "deep" copy onto TClonesArray (w/ placement new)
    AcdTkrIntersection* addAcdTkrIntersection(AcdTkrIntersection& inter);
    
    inline UInt_t nAcdIntersections() const { return m_acdTkrIntersectionCol != 0 ? 
						m_acdTkrIntersectionCol->GetEntriesFast() : 0; };
  
    inline const AcdTkrIntersection* getAcdTkrIntersection(UInt_t i) const { 
      return i < nAcdIntersections() ? 
	static_cast<const AcdTkrIntersection*>((*m_acdTkrIntersectionCol)[i]) : 0;  
    }  
    
private:

    /// Total energy in MeV deposited in the whole ACD system
    Double_t m_totEnergy;
    /// Distance of Closest Approach for the reconstructed gamma, 
    /// if there is one
    Double_t m_gammaDoca;
    /// Minimum Distance of Closest Approach for all tracks and all ACD tiles
    Double_t m_doca;
    /// New Bill Atwood DOCA calculation using edge of tiles
    Double_t m_actDist;
    /// Total number of ACD tiles above threshold
    Int_t m_tileCount;
    /// Collection of distance of closest approach calculations
    /// for each side row of the ACD
    vector<Double_t> m_rowDocaCol;
    /// Collection of Active Distance calc for each side row of ACD
    vector<Double_t> m_rowActDistCol;

    /// record of the tile with the minimum Distance of Closest Approach
    AcdId m_minDocaId;
    
    // Stores reconstructed energy per ACD digi
    vector<AcdId> m_idCol;
    vector<Double_t> m_energyCol;  

    /// record of the tile with the maximum Active Distance
    AcdId m_maxActDistId;
    /// Total MC Energy (MeV) deposited in the ribbons
    Double_t m_totRibbonEnergy;
    /// Total number of hit ribbons
    Int_t m_ribbonCount;
 
    Double32_t m_ribbonActDist;
    AcdId m_ribbonActDistId;

    // Store the expected track intersection points
    TClonesArray *m_acdTkrIntersectionCol;

    ClassDef(AcdRecon,8) // Acd Reconstruction data
};

#endif
