
#ifndef Root_CalRangeRec_H
#define Root_CalRangeRec_H 1

#include "TObject.h"
#include "TVector3.h"
#include "digiRootData/CalXtalId.h"

/** @class CalRangeRecData  
* @brief Position and reconstructed energies for both faces of Xtal for Cal                                
*
* @author  A.Chekhtman
*/
class CalRangeRecData : public TObject {  
    
public:
    CalRangeRecData() {};
    
    CalRangeRecData(Char_t rangeP, Double_t eneP, Char_t rangeM, Double_t eneM) :
    m_rangeP(rangeP), 
        m_eneP(eneP), 
        m_rangeM(rangeM), 
        m_eneM(eneM),
        m_pos(TVector3(0.,0.,0.))
    {};
    
    ~CalRangeRecData() {};
    
    void Clear(Option_t *option="");

    void Print(Option_t *option="");

    void initialize(const TVector3 &pos) {
        m_pos = pos;
    }
    
    TVector3 getPosition() const { return m_pos;}
    
    /// retrieve energy from specified face
    Double_t getEnergy(CalXtalId::XtalFace face) const;
    
    /// retrieve energy range from specified face
    Char_t getRange(CalXtalId::XtalFace face) const;
    
    
private:
    
    Double_t m_eneP, m_eneM;
    TVector3 m_pos;
    Char_t  m_rangeP, m_rangeM;
    
    ClassDef(CalRangeRecData,1)  // Reconstructed CAL xtals
};



#endif