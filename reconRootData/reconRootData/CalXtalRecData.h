
#ifndef Root_CalXtalRecCol_H
#define Root_CalXtalRecCol_H 1


// Include files
// Root Includes
#include "CalRangeRecData.h"
#include "TObject.h"
#include "TVector3.h"
#include <iostream>
#include <vector>
#include "digiRootData/CalXtalId.h"


/*!
//------------------------------------------------------------------------------
//
// \class   CalXtalRecData        
//  
// \brief reconstructed data for a calorimeter crystal                                
//              
// Author:  A.Chekhtman, Apr,17,2002
//
//------------------------------------------------------------------------------
*/

class CalXtalRecData :  public TObject { 
    
public:
    
    
    CalXtalRecData() {};
    
    //CalDigi(CalTrigMode mode, CalXtalId CalXtalId, ObjectVector<CalXtalReadout> readout) : 
    //    m_mode(mode),
    //        m_XtalId(CalXtalId),
    //        m_readout(readout)
    //{};
    
    /// Destructor
    virtual ~CalXtalRecData() { };
    
    /// Retrieve readout mode
     const CalXtalId::CalTrigMode getMode() const { return m_mode; };
    void setMode(CalXtalId::CalTrigMode m) { m_mode = m; };

    /// Retrieve Xtal identifier
     const CalXtalId getPackedId() const { return m_XtalId; };
     void setPackedId(CalXtalId id) { m_XtalId = id; };
    
     void addRangeRecData(CalRangeRecData r) { m_RecData.push_back(r); } ;
    
    /// Retrieve energy range for selected face and readout
     Char_t getRange(Short_t readoutIndex, CalXtalId::XtalFace face) const;
    
    /// Retrieve energy for selected face and readout
     Double_t getEnergy(Short_t readoutIndex, CalXtalId::XtalFace face) const;
    
    
    /// Retrieve average energy of two faces for the best range
     Double_t getEnergy();
    
    /// Retrieve the position for the best range
     TVector3 getPosition()
    {
        return getRangeRecData(0)->getPosition();
    }
    
    /// Retrieve reconstructed data from both ends of selected readout
    CalRangeRecData* getRangeRecData(short readoutIndex);
    
    /// Retrieve pulse height from selected range
    Short_t getEnergySelectedRange(char range, CalXtalId::XtalFace face) const;
    
    
    
private:
    
    /// Cal readout mode is based on trigger type
    CalXtalId::CalTrigMode m_mode;
    /// Cal ID
    CalXtalId m_XtalId;
    /// ranges and pulse heights
    std::vector<CalRangeRecData> m_RecData;
    
    ClassDef(CalXtalRecData,1)
};

#endif