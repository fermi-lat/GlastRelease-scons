#ifndef ROOT_ACDRECON_H
#define ROOT_ACDRECON_H

#include "TObject.h"
#include "digiRootData/AcdId.h"
#include <vector>
#ifdef WIN32
using namespace std;
#endif

/** 
* @class AcdRecon
*
* @brief Root container for the reconstructed ACD information 
* 
* 
* @author Heather Kelly
*
* $Header$
*/

class AcdRecon : public TObject
{
    
public:
    AcdRecon();
 
    AcdRecon(Double_t e, Int_t count, Double_t gDoca, Double_t doca, Double_t actDist,
            const AcdId &minDocaId, const std::vector<Double_t> &rowDoca)       
            : m_totEnergy(e),
            m_tileCount(count),
            m_gammaDoca(gDoca),
            m_doca(doca),
            m_actDist(actDist),
            m_minDocaId(minDocaId),
            m_rowDocaCol(rowDoca)
        {};
     
    virtual ~AcdRecon();
    
    void initialize(Double_t e, Int_t count, Double_t gDoca, Double_t doca, Double_t actDist,
            const AcdId &minDocaId, const std::vector<Double_t> &rowDoca);

    /// overload initialize for interactive root where we cannot use std::vector
    void initialize(Double_t e, Int_t count, Double_t gDoca, Double_t doca, Double_t actDist,
            const AcdId &minDocaId);
    
    void Clear(Option_t *option="");
    
    void Print(Option_t *option="") const;
    
    inline const Double_t getEnergy() const { return m_totEnergy; };
    inline const Int_t getTileCount() const { return m_tileCount; };
    inline const Double_t getGammaDoca() const { return m_gammaDoca; };
    inline const Double_t getDoca() const { return m_doca; };
    inline const Double_t getActiveDist() const { return m_actDist; };
    inline const AcdId& getMinDocaId() const { return m_minDocaId; };
    inline const std::vector<Double_t>& getRowDocaCol() const { return m_rowDocaCol; };
    inline const Double_t getRowDoca(UInt_t i) const { 
        return ((i < m_rowDocaCol.size()) ? m_rowDocaCol[i] : -1.); };
    inline void addRowDoca(Double_t val) { m_rowDocaCol.push_back(val); };
    //inline const std::map<AcdId, Double_t>& getEnergyCol() const { return m_energyCol; };
    
private:
    /// Total energy in MeV deposited in the whole ACD system
    Double_t m_totEnergy;
    /// Distance of Closest Approach for the reconstructed gamme, 
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
    std::vector<Double_t> m_rowDocaCol;
    
    // record of the tile with the minimum Distance of Closest Approach
    AcdId m_minDocaId;
    
    /// Stores reconstructed energy per ACD digi
    //std::map<AcdId, Double_t> m_energyCol;
    
    ClassDef(AcdRecon,1) // Acd Reconstruction data
};

#endif
