
// Include files
// Root Includes
#include "reconRootData/CalXtalRecData.h"

ClassImp(CalXtalRecData)

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

    /// Retrieve energy for selected face and readout
Double_t CalXtalRecData::getEnergy(Short_t readoutIndex, CalXtalId::XtalFace face) const
    {
        return (readoutIndex < m_RecData.size()) ? ((m_RecData[readoutIndex])).getEnergy(face) : (short)-1;
    };
    
    
    /// Retrieve average energy of two faces for the best range
 Double_t CalXtalRecData::getEnergy()
    {
        return (getEnergy(0,CalXtalId::POS)
            +getEnergy(0,CalXtalId::NEG))/2;
    };
    
    /// Retrieve energy range for selected face and readout
 Char_t CalXtalRecData::getRange(Short_t readoutIndex, CalXtalId::XtalFace face) const 
 {
        return (readoutIndex < m_RecData.size()) ? ((m_RecData[readoutIndex])).getRange(face) : (char)-1;
    }

 
 /// Retrieve reconstructed data from both ends of selected readout
 CalRangeRecData* CalXtalRecData::getRangeRecData(short readoutIndex)
    {
        //return ((readoutIndex < m_readout.size()) ? m_readout[readoutIndex] : 0);
        if ( readoutIndex < m_RecData.size() )
            return &(m_RecData[readoutIndex]);
        else
            return 0;
        
    };
    
    /// Retrieve pulse height from selected range
  short CalXtalRecData::getEnergySelectedRange(char range, CalXtalId::XtalFace face) const
    {
        char nRanges = (char)m_RecData.size();
        if (nRanges == 1)
            return (range == ((m_RecData[0])).getRange(face)) ? ((m_RecData[0])).getEnergy(face) : (short)-1;
        else
            return ((m_RecData[(nRanges + range - ((m_RecData[0])).getRange(face)) % nRanges])).getEnergy(face);
    };

