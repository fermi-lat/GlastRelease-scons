///////////////////////////////////////////////////////////////////////////
//                                                                       
// The L1T class contains the information about the Level 1 Trigger count
// for an event.  
//                                                                       
///////////////////////////////////////////////////////////////////////////

#include "digiRootData/L1T.h"
#include "TClass.h"

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

/// Implement the streamer ourselves for now...to take advantage
/// of schema evolution - and to allow our TBEvent class to handle
/// both old (<= Root v2.25) Root files, and new (>= Root v3.00) files
void L1T::Streamer(TBuffer &R__b)
{
   // Stream an object of class L1T.

   if (R__b.IsReading()) {
      UInt_t R__s, R__c;
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); 
      if (R__v > 1) 
      { 
          L1T::Class()->ReadBuffer(R__b, this, R__v, R__s, R__c);
          return;
      }
      /// Old Versions      
      TObject::Streamer(R__b);
      R__b >> m_trigCount;
      R__b >> m_trigTimeValue;
      R__b >> m_xCapture;
      R__b >> m_yCapture;
      R__b >> m_vetoCapture;
      R__b >> m_deadTimeReason;
      R__b >> m_deadTime;
   } else {
       L1T::Class()->WriteBuffer(R__b, this);
   }
}
