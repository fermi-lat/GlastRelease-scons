///////////////////////////////////////////////////////////////////////////
//									 
// The ESAPID class contains the information about particle identification
// and beam status.
//									 
///////////////////////////////////////////////////////////////////////////

#include "digiRootData/ESAPID.h"
//#include "TClass.h"

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

