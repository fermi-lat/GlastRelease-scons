
#ifndef GlastEvent_CalDigi_H
#define GlastEvent_CalDigi_H 1


// Include files
#include <iostream>
#include <vector>
#include "idents/CalXtalId.h"

#include "GaudiKernel/Kernel.h"
#include "GaudiKernel/StreamBuffer.h"
#include "GaudiKernel/ContainedObject.h"
#include "GaudiKernel/ObjectVector.h"
#include "GlastEvent/TopLevel/Definitions.h"





/*!
//------------------------------------------------------------------------------
//
// \class   CalDigi        
//  
// \brief Digitizations for Cal                                
//              
// Author:  J. Eric Grove, 23 Feb 2001
//
//------------------------------------------------------------------------------
*/

extern const CLID& CLID_CalDigi;

class CalDigi : virtual public ContainedObject { 
    
public:
    
/*!
//------------------------------------------------------------------------------
//
// \class   CalXtalReadout        
//  
// \brief   Pulse heights and energy range for both faces of Xtal for Cal                                
//              
// Author:  J. Eric Grove, 23 Feb 2001
//
//------------------------------------------------------------------------------
    */
    class CalXtalReadout {  // : virtual public ContainedObject  { 
        
    public:
        CalXtalReadout(char rangeP, short adcP, char rangeM, short adcM) :
          m_rangeP(rangeP), 
              m_adcP(adcP), 
              m_rangeM(rangeM), 
              m_adcM(adcM)
          {};
          
          /// Destructor
          ~CalXtalReadout() {};
          
          
          // retrieve pulse height from specified face
          inline short getAdc(idents::CalXtalId::XtalFace face) const {return face == idents::CalXtalId::POS ? m_adcP : m_adcM;};
          
          // retrieve energy range from specified face
          inline char getRange(idents::CalXtalId::XtalFace face) const {return face == idents::CalXtalId::POS ? m_rangeP : m_rangeM;};
          
          
    private:
        
        short m_adcP, m_adcM;
        char  m_rangeP, m_rangeM;
        
    };

    typedef std::vector<CalXtalReadout> CalXtalReadoutVector;
    
    /// shifts and masks for packed readout of energy range and Adc value
    enum {POS_OFFSET = 14,						// shift for POSitive face
        RANGE_OFFSET = 12, RANGE_MASK = 0x3000,		// energy range bits
        ADC_VAL_MASK = 0xfff};						// Adc value bits
    
    CalDigi() {};
    
    //CalDigi(CalTrigMode mode, idents::CalXtalId CalXtalId, ObjectVector<CalXtalReadout> readout) : 
    //    m_mode(mode),
    //        m_XtalId(CalXtalId),
    //        m_readout(readout)
    //{};
    
    /// Destructor
    virtual ~CalDigi() { };
    
    /// Retrieve readout mode
    inline const idents::CalXtalId::CalTrigMode getMode() const { return m_mode; }
    inline void setMode(idents::CalXtalId::CalTrigMode m) { m_mode = m; }
    
    /// Retrieve Xtal identifier
    inline const idents::CalXtalId getPackedId() const { return m_XtalId; }
    inline void setPackedId(idents::CalXtalId id) { m_XtalId = id; }
    
    inline void addReadout(CalXtalReadout r) { m_readout.push_back(r); } 
	const CalXtalReadoutVector& getReadoutVector() const { return m_readout;}	
    
    /// Retrieve energy range for selected face and readout
    inline char getRange(short readoutIndex, idents::CalXtalId::XtalFace face) const
    {
        return (readoutIndex < m_readout.size()) ? ((m_readout[readoutIndex])).getRange(face) : (char)-1;
    }
    
    /// Retrieve pulse height for selected face and readout
    inline short getAdc(short readoutIndex, idents::CalXtalId::XtalFace face) const
    {
        return (readoutIndex < m_readout.size()) ? ((m_readout[readoutIndex])).getAdc(face) : (short)-1;
    }
    
    /// Retrieve ranges and pulse heights from both ends of selected readout
    inline const CalXtalReadout* getXtalReadout(short readoutIndex)
    {
        //return ((readoutIndex < m_readout.size()) ? m_readout[readoutIndex] : 0);
        if ( readoutIndex < m_readout.size() )
            return &(m_readout[readoutIndex]);
        else
            return 0;
        
    }
    
    /// Retrieve pulse height from selected range
    inline short getAdcSelectedRange(char range, idents::CalXtalId::XtalFace face) const
    {
        char nRanges = (char)m_readout.size();
        if (nRanges == 1)
            return (range == ((m_readout[0])).getRange(face)) ? ((m_readout[0])).getAdc(face) : (short)-1;
        else
            return ((m_readout[(nRanges + range - ((m_readout[0])).getRange(face)) % nRanges])).getAdc(face);
    }

    

private:
    
    /// Cal readout mode is based on trigger type
    idents::CalXtalId::CalTrigMode m_mode;
    /// Cal ID
    idents::CalXtalId m_XtalId;
    /// ranges and pulse heights
    CalXtalReadoutVector m_readout;
    
};

typedef ObjectVector<CalDigi> CalDigiCol;


#endif