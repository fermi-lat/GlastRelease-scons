// File and Version Information:
//      $Header$
//
// Description:
//  Implement the ObfFilterStatus class for storing to root the output of the OnboardFilter
//                                                
// Author:  Tracy Usher                    

#include "digiRootData/ObfFilterStatus.h"
#include "Riostream.h"

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
    m_obfStatusCol.AddAt(0,HFCFilter);
    m_obfStatusCol.AddAt(0,MipFilter);
    m_obfStatusCol.AddAt(0,DFCFilter);

    return;
}


ObfFilterStatus& ObfFilterStatus::operator =(const ObfFilterStatus& rhs)
{
    m_obfStatusCol[GammaFilter] = rhs.m_obfStatusCol[GammaFilter];
    m_obfStatusCol[HFCFilter]   = rhs.m_obfStatusCol[HFCFilter];
    m_obfStatusCol[MipFilter]   = rhs.m_obfStatusCol[MipFilter];
    m_obfStatusCol[DFCFilter]   = rhs.m_obfStatusCol[DFCFilter];

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
    if (m_obfStatusCol[HFCFilter])   m_obfStatusCol[HFCFilter]->Print(option);
    if (m_obfStatusCol[MipFilter])   m_obfStatusCol[MipFilter]->Print(option);
    if (m_obfStatusCol[DFCFilter])   m_obfStatusCol[DFCFilter]->Print(option);
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

// This for backward compatibility...
UInt_t ObfGammaStatus::getStatusHi() const {return m_status >> enums::GFC_STATUS_V_TKR_LT_2_ELO;}
UInt_t ObfGammaStatus::getStatusLo() const {return m_status & 0xFFFF;}

// If msb is set below then event is to be vetoed
UInt_t ObfGammaStatus::getStatus32() const {return m_status;}
UInt_t ObfGammaStatus::getVetoMask() const {return enums::GFC_STATUS_M_VETOES;}
UInt_t ObfGammaStatus::getVetoBit()  const {return enums::GFC_STATUS_M_VETOED;}
// root Clear...
void ObfGammaStatus::Clear(Option_t *option)
{
    m_status = 0;
    return;
}

// root print...
void ObfGammaStatus::Print(Option_t* option) const 
{
    TObject::Print(option);

    std::cout << "Gamma Filter status:" << std::hex << m_status << std::dec << std::endl;
    return;
}

// Implement the ObfHFCStatus class here

#include "enums/XfcHFCStatus.h"
ClassImp(ObfHFCStatus)

// This for backward compatibility...
UInt_t ObfHFCStatus::getStatusHi() const {return m_status >> 16;}
UInt_t ObfHFCStatus::getStatusLo() const {return m_status & 0xFFFF;}

// If msb is set below then event is to be vetoed
UInt_t ObfHFCStatus::getStatus32() const {return m_status;}
UInt_t ObfHFCStatus::getVetoMask() const {return enums::HFC_STATUS_M_VETO_DEF;}
UInt_t ObfHFCStatus::getVetoBit()  const {return enums::HFC_STATUS_M_VETOED;}
// root Clear...
void ObfHFCStatus::Clear(Option_t *option)
{
    m_status = 0;
    return;
}

// root print...
void ObfHFCStatus::Print(Option_t* option) const 
{
    TObject::Print(option);

    std::cout << "HFC status:" << std::hex << m_status << std::dec << std::endl;
    return;
}
// Implement the ObfHFCStatus class here

#include "enums/XfcMipStatus.h"
ClassImp(ObfMipStatus)

// This for backward compatibility...
UInt_t ObfMipStatus::getStatusHi() const {return m_status >> 16;}
UInt_t ObfMipStatus::getStatusLo() const {return m_status & 0xFFFF;}

// If msb is set below then event is to be vetoed
UInt_t ObfMipStatus::getStatus32() const {return m_status;}
UInt_t ObfMipStatus::getVetoMask() const {return enums::MFC_STATUS_M_VETO_DEF;}
UInt_t ObfMipStatus::getVetoBit()  const {return enums::MFC_STATUS_M_VETOED;}
// root Clear...
void ObfMipStatus::Clear(Option_t *option)
{
    m_status = 0;
    return;
}

// root print...
void ObfMipStatus::Print(Option_t* option) const 
{
    TObject::Print(option);

    std::cout << "Mip Filter status:" << std::hex << m_status << std::dec << std::endl;
    return;
}


// Implement the ObfDFCStatus class here

#include "enums/XfcDFCStatus.h"
ClassImp(ObfDFCStatus)

// This for backward compatibility...
UInt_t ObfDFCStatus::getStatusHi() const {return m_status >> 16;}
UInt_t ObfDFCStatus::getStatusLo() const {return m_status & 0xFFFF;}

// If msb is set below then event is to be vetoed
UInt_t ObfDFCStatus::getStatus32() const {return m_status;}
UInt_t ObfDFCStatus::getVetoMask() const {return enums::DFC_STATUS_M_VETO_DEF;}
UInt_t ObfDFCStatus::getVetoBit()  const {return enums::DFC_STATUS_M_VETOED;}
// root Clear...
void ObfDFCStatus::Clear(Option_t *option)
{
    m_status = 0;
    return;
}

// root print...
void ObfDFCStatus::Print(Option_t* option) const 
{
    TObject::Print(option);

    std::cout << "DFC status:" << std::hex << m_status << std::dec << std::endl;
    return;
}
