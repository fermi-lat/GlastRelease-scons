
#ifndef CalRecon_CalXtalRecData_H
#define CalRecon_CalXtalRecData_H 1

#include <iostream>
#include <vector>
#include "idents/CalXtalId.h"

#include "GaudiKernel/Kernel.h"
#include "GaudiKernel/StreamBuffer.h"

#include "GaudiKernel/ObjectVector.h"
#include "GaudiKernel/ContainedObject.h"
#include "Event/TopLevel/Definitions.h"
#include "geometry/Point.h"



extern const CLID& CLID_CalXtalRecData;

namespace Event 
{

/** @class   CalXtalRecData        
 * @brief reconstructed data for a calorimeter crystal                                            
 * @author  A.Chekhtman
 * $Header$
*/
    class CalXtalRecData : virtual public ContainedObject { 
        
    public:
        
    /** @class   CalRangeRecdata        
     * @brief   position, reconstructed energies for both faces of Xtal for Cal            
     * @author  A.Chekhtman
     */
        class CalRangeRecData {  
            
        public:
            CalRangeRecData(char rangeP, double eneP, char rangeM, double eneM) :
              m_rangeP(rangeP), 
                  m_eneP(eneP), 
                  m_rangeM(rangeM), 
                  m_eneM(eneM),
                  m_pos(Point(0.,0.,0.))
              {};
              
              ~CalRangeRecData() {};
              
              void setPosition (Point pos) { m_pos = pos;}
              Point getPosition() const { return m_pos;}
              
              // retrieve energy from specified face
              inline double getEnergy(idents::CalXtalId::XtalFace face) const {return face == idents::CalXtalId::POS ? m_eneP : m_eneM;};
              
              // retrieve energy range from specified face
              inline char getRange(idents::CalXtalId::XtalFace face) const {return face == idents::CalXtalId::POS ? m_rangeP : m_rangeM;};
              
              
        private:
            
            double m_eneP, m_eneM;
            Point m_pos;
            char  m_rangeP, m_rangeM;
            
        };
        
        
        CalXtalRecData() {};
        
        CalXtalRecData(idents::CalXtalId::CalTrigMode mode, idents::CalXtalId CalXtalId) : 
        m_mode(mode), m_xtalId(CalXtalId){};
        
        virtual ~CalXtalRecData() { };
        
        void initialize (idents::CalXtalId::CalTrigMode m, idents::CalXtalId id)
        {m_mode = m; m_xtalId = id; }
        
        /// Retrieve readout mode
        inline const idents::CalXtalId::CalTrigMode getMode() const { return m_mode; };
        
        /// Retrieve Xtal identifier
        inline const idents::CalXtalId getPackedId() const { return m_xtalId; };
        
        inline void addRangeRecData(CalRangeRecData r) { m_RecData.push_back(r); } ;
        
        /// Retrieve energy range for selected face and readout
        inline char getRange(short readoutIndex, idents::CalXtalId::XtalFace face) const
        {
            return (readoutIndex < m_RecData.size()) ? ((m_RecData[readoutIndex])).getRange(face) : (char)-1;
        }
        
        /// Retrieve energy for selected face and readout
        inline double getEnergy(short readoutIndex, idents::CalXtalId::XtalFace face) const
        {
            return (readoutIndex < m_RecData.size()) ? ((m_RecData[readoutIndex])).getEnergy(face) : (short)-1;
        }
        
        
        /// Retrieve average energy of two faces for the best range
        inline double getEnergy()
        {
            return (getEnergy(0,idents::CalXtalId::POS)
                +getEnergy(0,idents::CalXtalId::NEG))/2;
        }
        
        /// Retrieve the position for the best range
        inline Point getPosition()
        {
            return getRangeRecData(0)->getPosition();
        }
        
        /// Retrieve reconstructed data from both ends of selected readout
        inline CalRangeRecData* getRangeRecData(short readoutIndex)
        {
            //return ((readoutIndex < m_readout.size()) ? m_readout[readoutIndex] : 0);
            if ( readoutIndex < m_RecData.size() )
                return &(m_RecData[readoutIndex]);
            else
                return 0;
            
        }
        
        /// Retrieve pulse height from selected range
        inline double getEnergySelectedRange(char range, idents::CalXtalId::XtalFace face) const
        {
            char nRanges = (char)m_RecData.size();
            if (nRanges == 1)
                return (range == ((m_RecData[0])).getRange(face)) ? ((m_RecData[0])).getEnergy(face) : -1.0;
            else
                return ((m_RecData[(nRanges + range - ((m_RecData[0])).getRange(face)) % nRanges])).getEnergy(face);
        }
        
        /// Serialize the object for writing
        //    virtual StreamBuffer& serialize( StreamBuffer& s ) const;
        /// Serialize the object for reading
        //    virtual StreamBuffer& serialize( StreamBuffer& s );
        /// Fill the ASCII output stream
        //    virtual std::ostream& fillStream( std::ostream& s ) const;
        
        
private:
    
    /// Cal readout mode is based on trigger type
    idents::CalXtalId::CalTrigMode m_mode;
    /// Cal ID
    idents::CalXtalId m_xtalId;
    /// ranges and pulse heights
    std::vector<CalRangeRecData> m_RecData;
    
};
typedef ObjectVector<CalXtalRecData> CalXtalRecCol;    

}


#endif
