// File and Version Information:
//      $Header$
//
// Description:
//  Implement the LpaFilterStatus class for storing to root the output of the LPA Event Handlers
//                                                
// Author:  Heather Kelly                   

#include "digiRootData/LpaHandler.h"
#include "Riostream.h"


ClassImp(ILpaHandler)


void ILpaHandler::initialize(UInt_t masterKey, UInt_t cfgKey, UInt_t cfgId, 
        enums::Lsf::RsdState state,
        enums::Lsf::LeakedPrescaler prescaler, UInt_t version, 
        enums::Lsf::HandlerId id, Bool_t has, UInt_t prescaleFactor) {
            m_masterKey = masterKey;
            m_cfgKey = cfgKey;
            m_cfgId = cfgId;
            m_state = state;
            m_prescaler = prescaler;
            m_version = version;
            m_id = id;
            m_has = has;
            setPrescaleFactor(prescaleFactor);
        }

 void ILpaHandler::Clear(Option_t* /*option*/) {
     m_type = enums::Lsf::Unknown;
     m_masterKey = 0xFFFFFFFF;
     m_cfgKey = 0xFFFFFFFF; 
     m_cfgId = 0xFFFFFFFF;
     m_state = enums::Lsf::INVALID;
     m_prescaler = enums::Lsf::UNSUPPORTED;
     m_version = 0;
     m_id = enums::Lsf::MaxHandlerIds;
     m_has = kFALSE;
     m_prescaleFactor = 0;
 }

 void ILpaHandler::Print(Option_t* option) const 
{
    TObject::Print(option);

    std::cout << "ILpaHandler:" << std::endl;
    return;
}


ClassImp(LpaHandler)

LpaHandler::LpaHandler() : m_gamma(0), m_hip(0), m_mip(0), m_dgn(0), m_pass(0)
   , m_asc(0)
{
    initialize();
}

LpaHandler::~LpaHandler() 
{    
	Clear();
}

void LpaHandler::initialize() 
{
    return;
}


LpaHandler& LpaHandler::operator =(const LpaHandler& rhs)
{
    if (rhs.m_gamma) m_gamma = new LpaGammaFilter(*(rhs.getGammaFilter()));
    if (rhs.m_hip) m_hip = new LpaHipFilter(*(rhs.getHipFilter()));
    if (rhs.m_mip) m_mip = new LpaMipFilter(*(rhs.getMipFilter()));
    if (rhs.m_dgn) m_dgn = new LpaDgnFilter(*(rhs.getDgnFilter()));
    if (rhs.m_pass) m_pass = new LpaPassthruFilter(*(rhs.getPassthruFilter()));
    return *this;
}

void LpaHandler::Clear(Option_t *) 
{
    if(m_gamma) { 
       delete m_gamma;
       m_gamma = 0;
    }
    if (m_hip) {
        delete m_hip;
        m_hip = 0;
    }
   if (m_mip) {
        delete m_mip;
        m_mip =0;
   }
   if (m_dgn) {
       delete m_dgn;
       m_dgn = 0;
   }
   if (m_pass) {
       delete m_pass;
       m_pass = 0;
   }
   if (m_asc) {
       delete m_asc;
       m_asc = 0;
   }
    return;
}

// dummy data, just for tests
void LpaHandler::Fake( Int_t /*ievent*/, UInt_t /*rank*/, Float_t /*randNum*/ ) {

    Clear() ;
}

//#define COMPARE_IN_RANGE(att) rootdatautil::CompareInRange(get ## att(),ref.get ## att(),#att)

Bool_t LpaHandler::CompareInRange( const LpaHandler & /*ref*/, const std::string & /*name*/ ) const {

    Bool_t result = true ;
   
    return result ;

}

void LpaHandler::Print(Option_t* /*option*/) const 
{
}

// Add results method
void LpaHandler::addHandler(const enums::Lsf::HandlerId &id, TObject* status)
{
    //m_lpaHandlerCol.AddAtAndExpand(status, id);
    if (id == enums::Lsf::GAMMA) m_gamma = dynamic_cast<LpaGammaFilter*>(status);
    else if (id == enums::Lsf::HIP) m_hip = dynamic_cast<LpaHipFilter*>(status);
    else if (id == enums::Lsf::MIP) m_mip = dynamic_cast<LpaMipFilter*>(status);
    else if (id == enums::Lsf::PASS_THRU) m_pass = dynamic_cast<LpaPassthruFilter*>(status);
    else if (id == enums::Lsf::DGN) m_dgn = dynamic_cast<LpaDgnFilter*>(status);
    return;
}

// Return results method
const ILpaHandler* LpaHandler::getHandler(const enums::Lsf::HandlerId &id) const
{
    if (id == enums::Lsf::GAMMA) return m_gamma;
    else if (id == enums::Lsf::HIP) return m_hip;
    else if (id == enums::Lsf::MIP) return m_mip;
    else if (id == enums::Lsf::PASS_THRU) return m_pass;
    else if (id == enums::Lsf::DGN) return m_dgn;
    else if (id == enums::Lsf::MaxHandlerIds) return 0;
    else return 0;
}

Int_t ILpaHandler::prescalerIndex() const {
  switch (m_prescaler){
  case enums::Lsf::UNSUPPORTED:
    // We really don't know what happened, maybe the state will tell us something
    break;
  case enums::Lsf::INPUT:
    // filter never saw the event, set a value larger than any specific prescaler
    return 33;
  case enums::Lsf::OUTPUT:
    // output prescale expired, this means the line prescalers are irrelevent
    // set a value larger that the first line
    return 32;      
  default:
    // one of the line prescales is asserted, return that 
    return m_prescaler;
  }

  switch (m_state) {
  case enums::Lsf::INVALID:  
    // We really don't know what happened, set a large value so this 
    // doesn't get into any of the other samples
    return 100;
  case enums::Lsf::IGNORED:
    // Filter never saw the event, set a value larger than any specific prescaler
    // should have been caught in the first switch above
    return 33;
  case enums::Lsf::PASSED:
    // If the event made it here it means that none of the bits were asserted
    // return a value smaller that the lowest bit
    return -1;
  case enums::Lsf::SUPPRESSED:
  case enums::Lsf::VETOED:
  case enums::Lsf::LEAKED:    
    // In all of these cases the event should have been caught by the first switch
    // but return the relevent prescaler, just in case
    break;
  }
  return m_prescaler;
}


// Implement the LpaGammaFilter class here

ClassImp(LpaGammaFilter)

// root Clear...
void LpaGammaFilter::Clear(Option_t* /*option*/)
{
    m_status = 0;
    return;
}

UInt_t LpaGammaFilter::getAllVetoBits() const { return (m_status & enums::GFC_STATUS_M_VETOES); }

// root print...
void LpaGammaFilter::Print(Option_t* option) const 
{
    TObject::Print(option);

    std::cout << "GAMMA Filter status:" << std::hex << m_status << std::dec << std::endl;
    return;
}

ClassImp(LpaGammaFilterV0)

ClassImp(LpaGammaFilterV1)

ClassImp(LpaGammaFilterV2)

ClassImp(LpaGammaFilterV3)

// Implement the LpaHipFilter class here

ClassImp(LpaHipFilter)

// root Clear...
void LpaHipFilter::Clear(Option_t* /*option*/)
{
    m_status = 0;
    return;
}

UInt_t LpaHipFilter::getAllVetoBits() const { return (m_status & enums::HFC_STATUS_M_VETO_DEF); }


// root print...
void LpaHipFilter::Print(Option_t* option) const 
{
    TObject::Print(option);

    std::cout << "HIP status:" << std::hex << m_status << std::dec << std::endl;
    return;
}

ClassImp(LpaHipFilterV0)

// Implement the LpaMipFilter class here

ClassImp(LpaMipFilter)


// root Clear...
void LpaMipFilter::Clear(Option_t* /*option*/)
{
    m_status = 0;
    return;
}

UInt_t LpaMipFilter::getAllVetoBits() const {return (m_status & enums::MFC_STATUS_M_VETO_DEF);}

// root print...
void LpaMipFilter::Print(Option_t* option) const 
{
    TObject::Print(option);

    std::cout << "MIP Filter status:" << std::hex << m_status << std::dec << std::endl;
    return;
}

ClassImp(LpaMipFilterV0)

// Implement the LpaDgnFilter class here

ClassImp(LpaDgnFilter)

// root Clear...
void LpaDgnFilter::Clear(Option_t* /*option*/)
{
    m_status = 0;
    return;
}

UInt_t LpaDgnFilter::getAllVetoBits() const {return (m_status & enums::DFC_STATUS_M_VETO_DEF);}


// root print...
void LpaDgnFilter::Print(Option_t* option) const 
{
    TObject::Print(option);

    std::cout << "DGN status:" << std::hex << m_status << std::dec << std::endl;
    return;
}

ClassImp(LpaDgnFilterV0)

ClassImp(LpaPassthruFilter)

// root Clear...
void LpaPassthruFilter::Clear(Option_t* /*option*/)
{
    m_status = 0;
    return;
}

// root print...
void LpaPassthruFilter::Print(Option_t* option) const 
{
    TObject::Print(option);

    std::cout << "Passthru status:" << std::hex << m_status << std::dec << std::endl;
    return;
}

ClassImp(LpaPassthruFilterV0)
