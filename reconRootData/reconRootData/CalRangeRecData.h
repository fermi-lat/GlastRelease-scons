
#ifndef Root_CalRangeRec_H
#define Root_CalRangeRec_H 1


// Include files
// Root Includes
#include "TObject.h"
#include "TVector3.h"
#include "digiRootData/CalXtalId.h"

/*!
//------------------------------------------------------------------------------
//
// \class   CalRangeRecdata        
//  
// \brief   position, reconstructed energies for both faces of Xtal for Cal                                
//              
// Author:  A.Chekhtman, Apr, 17, 2002
//
//------------------------------------------------------------------------------
*/
class CalRangeRecData : public TObject {  
    
public:
    CalRangeRecData() {};

    CalRangeRecData(char rangeP, Double_t eneP, char rangeM, Double_t eneM) :
      m_rangeP(rangeP), 
          m_eneP(eneP), 
          m_rangeM(rangeM), 
          m_eneM(eneM),
          m_pos(TVector3(0.,0.,0.))
      {};
      
      /// Destructor
      ~CalRangeRecData() {};
      
      void setPosition (TVector3 pos) { m_pos = pos;}
      TVector3 getPosition() const { return m_pos;}
      
      // retrieve energy from specified face
       Double_t getEnergy(CalXtalId::XtalFace face) const;
      
      // retrieve energy range from specified face
       Char_t getRange(CalXtalId::XtalFace face) const;
      
      
private:
    
    Double_t m_eneP, m_eneM;
    TVector3 m_pos;
    Char_t  m_rangeP, m_rangeM;
    ClassDef(CalRangeRecData,1)  
};



#endif