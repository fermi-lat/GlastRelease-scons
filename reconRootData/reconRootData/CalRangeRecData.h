
#ifndef Root_CalRangeRec_H
#define Root_CalRangeRec_H 1

#include "TObject.h"
#include "TVector3.h"
#include "commonRootData/idents/CalXtalId.h"

/** @class CalRangeRecData  
* @brief Position and reconstructed energies for both faces of Xtal for Cal                                
*
* @author  A.Chekhtman
*
* $Header$
*/
class CalRangeRecData : public TObject {  
    
public:
    CalRangeRecData() {};
    
    CalRangeRecData(Char_t rangeP, Double_t eneP, Char_t rangeM, Double_t eneM) :
    m_eneP(eneP), 
    m_eneM(eneM),
    m_pos(TVector3(0.,0.,0.)),
    m_rangeP(rangeP), 
    m_rangeM(rangeM)
    {};
    
    ~CalRangeRecData() {};
    
    void Clear(Option_t *option="");
    void Fake( Int_t ievent, UInt_t rank, Float_t randNum ) ; // for tests
    Bool_t CompareInRange( const CalRangeRecData &, const std::string & name = "" ) const ; // for tests
    void Print(Option_t *option="") const;

    void initialize(const TVector3 &pos) {
        m_pos = pos;
    }
    
    const TVector3& getPosition() const { return m_pos;}
    
    /// retrieve energy from specified face
    Double_t getEnergy(CalXtalId::XtalFace face) const;
    
    /// retrieve energy range from specified face
    Char_t getRange(CalXtalId::XtalFace face) const;
    
    
private:
    
        /// Energy associated with the plus end of the crystal
    Double32_t m_eneP;
        /// Energy associated with the minus end of the crystal
    Double32_t m_eneM;
        /// Position of the crystal
    TVector3 m_pos;
        /// Range associated with the plus end of the crystal
    Char_t  m_rangeP;
        /// Range associated with the minus end of the crystal
    Char_t m_rangeM;
    
    ClassDef(CalRangeRecData,2)  // Reconstructed CAL xtals
};



#endif
