///////////////////////////////////////////////////////////////////////////
//                                                                       
// The L1T class contains the information about the Level 1 Trigger count
// for an event.  
//                                                                       
///////////////////////////////////////////////////////////////////////////

#include "digiRootData/L1T.h"
//#include "TClass.h"

ClassImp(L1T)

///________________________________________________________________________
L1T::L1T( Int_t trigCount, Int_t trigTimeValue, Int_t xCapture, Int_t yCapture,
		 UShort_t vetoCapture, UChar_t deadTimeReason, UShort_t deadTime) :
	m_trigCount(trigCount), m_trigTimeValue(trigTimeValue), 
	m_xCapture(xCapture), m_yCapture(yCapture), m_vetoCapture(vetoCapture),
	m_deadTimeReason(deadTimeReason), m_deadTime(deadTime) {
}
//_________________________________________________________________________
L1T::L1T(){
  // Default constructor
}
//_________________________________________________________________________
L1T::~L1T(){
  // Destructor
}

