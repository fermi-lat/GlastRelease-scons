///////////////////////////////////////////////////////////////////////////
//                                                                       
// The L1T class contains the information about the Level 1 Trigger count
// for an event.  
//                                                                       
///////////////////////////////////////////////////////////////////////////

#include "digiRootData/L1T.h"

ClassImp(L1T)

///________________________________________________________________________
L1T::L1T( Int_t trigCount, Int_t trigTimeValue, Int_t xCapture, Int_t yCapture,
		 UShort_t vetoCapture, UChar_t deadTimeCause, UShort_t deadTime) :
	m_trigCount(trigCount), m_trigTimeValue(trigTimeValue), 
	m_xCapture(xCapture), m_yCapture(yCapture), m_vetoCapture(vetoCapture),
	m_deadTimeCause(deadTimeCause), m_deadTime(deadTime) {
}
//_________________________________________________________________________
L1T::L1T() : m_trigCount(0), m_trigTimeValue(0), m_xCapture(0), m_yCapture(0),
m_vetoCapture(0), m_deadTimeCause(0), m_deadTime(0)
{
  // Default constructor
}
//_________________________________________________________________________
L1T::~L1T(){
  // Destructor
}

void L1T::Clean() {
    m_trigCount = 0;  
    m_trigTimeValue = 0;
    m_xCapture = 0;
    m_yCapture = 0; 
    m_vetoCapture = 0;  
    m_deadTimeCause = 0; 
    m_deadTime = 0;     
}
