// File and Version Information:
//      $Header$
//
// Description:
//  Implement the ObfFilterStatus class for storing to root the output of the OnboardFilter
//                                                
// Author:  Tracy Usher                    

#include "digiRootData/ObfFilterStatus.h"
#include "Riostream.h"

// Use this to define the "state" of each filter
#include "enums/Lsf.h"
#include <map>

ClassImp(ObfFilterStatus)

ObfFilterStatus::ObfFilterStatus() 
{
    initialize();
}

ObfFilterStatus::~ObfFilterStatus() 
{    
	Clear();
}

void ObfFilterStatus::initialize() 
{
    m_obfStatusCol.Clear();
    m_obfStatusCol.AddAt(0,GammaFilter);
    m_obfStatusCol.AddAt(0,HipFilter);
    m_obfStatusCol.AddAt(0,MipFilter);
    m_obfStatusCol.AddAt(0,DgnFilter);

    return;
}


ObfFilterStatus& ObfFilterStatus::operator =(const ObfFilterStatus& rhs)
{
    m_obfStatusCol[GammaFilter] = rhs.m_obfStatusCol[GammaFilter];
    m_obfStatusCol[HipFilter]   = rhs.m_obfStatusCol[HipFilter];
    m_obfStatusCol[MipFilter]   = rhs.m_obfStatusCol[MipFilter];
    m_obfStatusCol[DgnFilter]   = rhs.m_obfStatusCol[DgnFilter];

    m_obfStatusCol.SetOwner(kTRUE);

    return *this;
}

void ObfFilterStatus::Clear(Option_t *) 
{
    initialize();
    return;
}

// dummy data, just for tests
void ObfFilterStatus::Fake( Int_t /*ievent*/, UInt_t rank, Float_t randNum ) {

    Clear() ;
/*    
    Float_t f = Float_t(rank) ;
    Float_t fr = f*randNum ;

    ObfFilterStatusPoint* point = new ObfFilterStatusPoint();
    point->Fake(0, rank, randNum);
    addMcPoint(point);

    point = new ObfFilterStatusPoint();
    point->Fake(0, 2*rank, randNum);
    addMcPoint(point);
*/
}

//#define COMPARE_IN_RANGE(att) rootdatautil::CompareInRange(get ## att(),ref.get ## att(),#att)

Bool_t ObfFilterStatus::CompareInRange( const ObfFilterStatus & ref, const std::string & name ) const {

    Bool_t result = true ;
    
//    result = COMPARE_IN_RANGE(ObfFilterStatusCharge) && result 

//    const McParticle *myPart = getMcParticle();
//    const McParticle *refPart = ref.getMcParticle();
//    result = rootdatautil::CompareInRange(*myPart,*refPart,"McParticle") && result ;


    //result = COMPARE_IN_RANGE(McParticle) && result ;  need to pass reference
/*
    result = rootdatautil::TObjArrayCompareInRange<ObfFilterStatusPoint>(getMcPointCol(),ref.getMcPointCol(),"McPoint") && result ;
      
    if (!result) {
        if ( name == "" ) {
            std::cout<<"Comparison ERROR for "<<ClassName()<<std::endl ;
        }
        else {
            std::cout<<"Comparison ERROR for "<<name<<std::endl ;
        }
    }
*/
    return result ;

}

void ObfFilterStatus::Print(Option_t* option) const 
{
    if (m_obfStatusCol[GammaFilter]) m_obfStatusCol[GammaFilter]->Print(option);
    if (m_obfStatusCol[HipFilter])   m_obfStatusCol[HipFilter]->Print(option);
    if (m_obfStatusCol[MipFilter])   m_obfStatusCol[MipFilter]->Print(option);
    if (m_obfStatusCol[DgnFilter])   m_obfStatusCol[DgnFilter]->Print(option);
}

// Add results method
//void ObfFilterStatus::addFilterStatus(FilterKeys key, IObfStatus* status)
void ObfFilterStatus::addFilterStatus(FilterKeys key, TObject* status)
{
    m_obfStatusCol[key] = status;

    return;
}

// Return results method
const IObfStatus* ObfFilterStatus::getFilterStatus(FilterKeys key) const
{
    IObfStatus* status  = 0;
    TObject*    statObj = m_obfStatusCol[key];
    if (m_obfStatusCol[key]) status = dynamic_cast<IObfStatus*>(m_obfStatusCol[key]);

    return status;
}

// Implement the ObfGammaStatus class here

#include "enums/EfcGammaStatus.h"
ClassImp(ObfGammaStatus)

// If msb is set below then event is to be vetoed
UInt_t ObfGammaStatus::getStatusWord() const {return m_status;}
UInt_t ObfGammaStatus::getVetoMask() const {return enums::GFC_STATUS_M_VETOES;}
UInt_t ObfGammaStatus::getVetoBit()  const {return enums::GFC_STATUS_M_VETOED;}
// root Clear...
void ObfGammaStatus::Clear(Option_t *option)
{
    m_status = 0;
    return;
}

UInt_t ObfGammaStatus::getEnergyInLeus() const
{
     UInt_t gamEnergy = m_energy & 0x00FFFFFF;

     // drag the sign bit up and then down...
     gamEnergy <<= 8;
     gamEnergy >>= 8;

     return gamEnergy;
}

UInt_t ObfGammaStatus::getState() const
{
    std::map<unsigned int, int> bitToEnumMap;

    bitToEnumMap[0x00] = enums::Lsf::PASSED;
    bitToEnumMap[0x20] = enums::Lsf::SUPPRESSED;
    bitToEnumMap[0x40] = enums::Lsf::VETOED;
    bitToEnumMap[0x60] = enums::Lsf::LEAKED;

    if (bitToEnumMap.find(m_sb) != bitToEnumMap.end()) return bitToEnumMap[m_sb];
    else return enums::Lsf::INVALID;
}

// root print...
void ObfGammaStatus::Print(Option_t* option) const 
{
    TObject::Print(option);

    std::cout << "Gamma Filter status:" << std::hex << m_status << std::dec << std::endl;
    return;
}

// Implement the ObfHipStatus class here

#include "enums/XfcHFCStatus.h"
ClassImp(ObfHipStatus)

// If msb is set below then event is to be vetoed
UInt_t ObfHipStatus::getStatusWord() const {return m_status;}
UInt_t ObfHipStatus::getVetoMask() const {return enums::HFC_STATUS_M_VETO_DEF;}
UInt_t ObfHipStatus::getVetoBit()  const {return enums::HFC_STATUS_M_VETOED;}
// root Clear...
void ObfHipStatus::Clear(Option_t *option)
{
    m_status = 0;
    return;
}

UInt_t ObfHipStatus::getState() const
{
    std::map<unsigned int, int> bitToEnumMap;

    bitToEnumMap[0x00] = enums::Lsf::PASSED;
    bitToEnumMap[0x20] = enums::Lsf::SUPPRESSED;
    bitToEnumMap[0x40] = enums::Lsf::VETOED;
    bitToEnumMap[0x60] = enums::Lsf::LEAKED;

    if (bitToEnumMap.find(m_sb) != bitToEnumMap.end()) return bitToEnumMap[m_sb];
    else return enums::Lsf::INVALID;
}

// root print...
void ObfHipStatus::Print(Option_t* option) const 
{
    TObject::Print(option);

    std::cout << "HFC status:" << std::hex << m_status << std::dec << std::endl;
    return;
}

// Implement the ObfHipStatus class here
ClassImp(ObfHFCStatus)

// If msb is set below then event is to be vetoed
UInt_t ObfHFCStatus::getStatusWord() const {return m_status;}
UInt_t ObfHFCStatus::getVetoMask() const {return enums::HFC_STATUS_M_VETO_DEF;}
UInt_t ObfHFCStatus::getVetoBit()  const {return enums::HFC_STATUS_M_VETOED;}
// root Clear...
void ObfHFCStatus::Clear(Option_t *option)
{
    m_status = 0;
    return;
}

UInt_t ObfHFCStatus::getState() const
{
    std::map<unsigned int, int> bitToEnumMap;

    bitToEnumMap[0x00] = enums::Lsf::PASSED;
    bitToEnumMap[0x20] = enums::Lsf::SUPPRESSED;
    bitToEnumMap[0x40] = enums::Lsf::VETOED;
    bitToEnumMap[0x60] = enums::Lsf::LEAKED;

    if (bitToEnumMap.find(m_sb) != bitToEnumMap.end()) return bitToEnumMap[m_sb];
    else return enums::Lsf::INVALID;
}

// root print...
void ObfHFCStatus::Print(Option_t* option) const 
{
    TObject::Print(option);

    std::cout << "HFC status:" << std::hex << m_status << std::dec << std::endl;
    return;
}
// Implement the ObfHipStatus class here


#include "enums/XfcMipStatus.h"
ClassImp(ObfMipStatus)

// If msb is set below then event is to be vetoed
UInt_t ObfMipStatus::getStatusWord() const {return m_status;}
UInt_t ObfMipStatus::getVetoMask() const {return enums::MFC_STATUS_M_VETO_DEF;}
UInt_t ObfMipStatus::getVetoBit()  const {return enums::MFC_STATUS_M_VETOED;}
// root Clear...
void ObfMipStatus::Clear(Option_t *option)
{
    m_status = 0;
    return;
}

UInt_t ObfMipStatus::getState() const
{
    std::map<unsigned int, int> bitToEnumMap;

    bitToEnumMap[0x00] = enums::Lsf::PASSED;
    bitToEnumMap[0x20] = enums::Lsf::SUPPRESSED;
    bitToEnumMap[0x40] = enums::Lsf::VETOED;
    bitToEnumMap[0x60] = enums::Lsf::LEAKED;

    if (bitToEnumMap.find(m_sb) != bitToEnumMap.end()) return bitToEnumMap[m_sb];
    else return enums::Lsf::INVALID;
}

// root print...
void ObfMipStatus::Print(Option_t* option) const 
{
    TObject::Print(option);

    std::cout << "Mip Filter status:" << std::hex << m_status << std::dec << std::endl;
    return;
}


// Implement the ObfDgnStatus class here

#include "enums/XfcDFCStatus.h"
ClassImp(ObfDgnStatus)

// If msb is set below then event is to be vetoed
UInt_t ObfDgnStatus::getStatusWord() const {return m_status;}
UInt_t ObfDgnStatus::getVetoMask() const {return enums::DFC_STATUS_M_VETO_DEF;}
UInt_t ObfDgnStatus::getVetoBit()  const {return enums::DFC_STATUS_M_VETOED;}
// root Clear...
void ObfDgnStatus::Clear(Option_t *option)
{
    m_status = 0;
    return;
}

UInt_t ObfDgnStatus::getState() const
{
    std::map<unsigned int, int> bitToEnumMap;

    bitToEnumMap[0x00] = enums::Lsf::PASSED;
    bitToEnumMap[0x20] = enums::Lsf::SUPPRESSED;
    bitToEnumMap[0x40] = enums::Lsf::VETOED;
    bitToEnumMap[0x60] = enums::Lsf::LEAKED;

    if (bitToEnumMap.find(m_sb) != bitToEnumMap.end()) return bitToEnumMap[m_sb];
    else return enums::Lsf::INVALID;
}

// root print...
void ObfDgnStatus::Print(Option_t* option) const 
{
    TObject::Print(option);

    std::cout << "DFC status:" << std::hex << m_status << std::dec << std::endl;
    return;
}

ClassImp(ObfDFCStatus)

// If msb is set below then event is to be vetoed
UInt_t ObfDFCStatus::getStatusWord() const {return m_status;}
UInt_t ObfDFCStatus::getVetoMask() const {return enums::DFC_STATUS_M_VETO_DEF;}
UInt_t ObfDFCStatus::getVetoBit()  const {return enums::DFC_STATUS_M_VETOED;}
// root Clear...
void ObfDFCStatus::Clear(Option_t *option)
{
    m_status = 0;
    return;
}

UInt_t ObfDFCStatus::getState() const
{
    std::map<unsigned int, int> bitToEnumMap;

    bitToEnumMap[0x00] = enums::Lsf::PASSED;
    bitToEnumMap[0x20] = enums::Lsf::SUPPRESSED;
    bitToEnumMap[0x40] = enums::Lsf::VETOED;
    bitToEnumMap[0x60] = enums::Lsf::LEAKED;

    if (bitToEnumMap.find(m_sb) != bitToEnumMap.end()) return bitToEnumMap[m_sb];
    else return enums::Lsf::INVALID;
}

// root print...
void ObfDFCStatus::Print(Option_t* option) const 
{
    TObject::Print(option);

    std::cout << "DFC status:" << std::hex << m_status << std::dec << std::endl;
    return;
}

