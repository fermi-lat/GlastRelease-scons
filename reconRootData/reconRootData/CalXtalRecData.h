
#ifndef Root_CalXtalRecCol_H
#define Root_CalXtalRecCol_H 1


#include "CalRangeRecData.h"
#include "TObject.h"
#include "TVector3.h"
#include <vector>
#ifdef WIN32
using namespace std;
#endif

#include "digiRootData/CalXtalId.h"


/** @class   CalXtalRecData        
* @brief Reconstructed data for one calorimeter crystal                                
*             
* @author A.Chekhtman
* $Header$
*/

class CalXtalRecData :  public TObject { 
    
public:
    
    
    CalXtalRecData() {};
        
    virtual ~CalXtalRecData() { };
    
    void initialize(CalXtalId::CalTrigMode m, CalXtalId id);

    void Clear(Option_t *option="");

    void Print(Option_t *option="") const;

    /// Retrieve readout mode
    const CalXtalId::CalTrigMode getMode() const { return m_mode; };
    
    /// Retrieve Xtal identifier
    const CalXtalId getPackedId() const { return m_xtalId; };
    
    /// Add a new element to the CalRangeRecData collection
    void addRangeRecData(CalRangeRecData r) { m_recData.push_back(r); } ;
    
    /// Retrieve energy range for selected face and readout
    Char_t getRange(Short_t readoutIndex, CalXtalId::XtalFace face) const;
    
    /// Retrieve energy for selected face and readout
    Double_t getEnergy(Short_t readoutIndex, CalXtalId::XtalFace face) const;
    
    
    /// Retrieve average energy of two faces for the best range
    Double_t getEnergy() const;
    
    /// Retrieve the position for the best range
    const TVector3& getPosition() const
    {
        return getRangeRecData(0)->getPosition();
    }
    
    /// Retrieve reconstructed data from both ends of selected readout
    const CalRangeRecData* getRangeRecData(Short_t readoutIndex) const;
    
    /// Retrieve pulse height from selected range
    Double_t getEnergySelectedRange(CalXtalId::AdcRange range, CalXtalId::XtalFace face) const;
    
    
    
private:
    
    /// Cal readout mode is based on trigger type
    CalXtalId::CalTrigMode m_mode;
    /// Cal ID
    CalXtalId m_xtalId;
    /// ranges and pulse heights
    std::vector<CalRangeRecData> m_recData;
    
    ClassDef(CalXtalRecData,1)
};

#endif