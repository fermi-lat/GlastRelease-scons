///////////////////////////////////////////////////////////////////////////
//                                                                       
// The ESAPID class contains the information about particle identification
// and beam status.
//                                                                       
///////////////////////////////////////////////////////////////////////////

#include "digiRootData/ESAPID.h"
#include "TClass.h"

ClassImp(ESAPID)

///________________________________________________________________________
ESAPID::ESAPID(){
  // Default constructor
    Clean();
}
//_________________________________________________________________________
ESAPID::~ESAPID(){
  // Destructor
}
//_________________________________________________________________________
void ESAPID::Clean() {
    m_Clean = kFALSE;
    m_XBeam = 0;
    m_YBeam = 0;
    m_NPart = 0;
    int i;
    for (i=0; i<10; i++) {
        m_ID[i] = 0;
        m_Mom[i] = 0;
    }
    
    for (i=0; i<5; i++) {
      m_BeamPb[i] = 0;
    }

    m_CADC = 0;

    for (i=0; i<4; i++) {
      m_TOF[i] = 0;
    }

    m_TOFC = 0;
    m_OneGoodP = kFALSE;
    m_OneGoodG = kFALSE;
    m_OneGoodE = kFALSE;
    m_OneGoodK = kFALSE;
    m_OneGoodPi = kFALSE;
    m_OneGoodH = kFALSE;
}

/// Implement the streamer ourselves for now...to take advantage
/// of schema evolution - and to allow our TBEvent class to handle
/// both old (<= Root v2.25) Root files, and new (>= Root v3.00) files
void ESAPID::Streamer(TBuffer &R__b)
{
   // Stream an object of class ESAPID.
   if (R__b.IsReading()) {
      UInt_t R__s, R__c;
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); 
      if (R__v > 1) 
      { 
          ESAPID::Class()->ReadBuffer(R__b, this, R__v, R__s, R__c);
          return;
      }
      TObject::Streamer(R__b);
      R__b >> m_Clean;
      R__b >> m_XBeam;
      R__b >> m_YBeam;
      R__b >> m_NPart;
      R__b.ReadStaticArray(m_ID);
      R__b.ReadStaticArray(m_Mom);
      R__b.ReadStaticArray(m_BeamPb);
      R__b >> m_CADC;
      R__b.ReadStaticArray(m_TOF);
      R__b >> m_TOFC;
      R__b >> m_OneGoodP;
      R__b >> m_OneGoodG;
      R__b >> m_OneGoodE;
      R__b >> m_OneGoodK;
      R__b >> m_OneGoodPi;
      R__b >> m_OneGoodH;
   } else {
       ESAPID::Class()->WriteBuffer(R__b, this);
   }
}