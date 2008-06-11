#ifndef obfFilterStatusRoot_H
#define obfFilterStatusRoot_H

#include "TObject.h"
#include "TObjArray.h"

/** @class IObfStatus
* @brief Interface/base class which describes access to the results of the various 
*        onboard filters. Currently, all filters output a 32 bit status word but then
*        can also include some number of other items. This interface class deals 
*        entirely with the status bits
*
*/
class IObfStatus
{
public:
    virtual ~IObfStatus() {}

    /// TObject functions
    virtual void Clear(Option_t *option ="") = 0;
    virtual void Print(Option_t *option="") const = 0;

    /// OBF functions
    virtual UInt_t  getStatusWord() const = 0;
    virtual UInt_t  getVetoMask()   const = 0;
    virtual UInt_t  getVetoBit()    const = 0;

    virtual UChar_t getFilterId()   const = 0;
    virtual UChar_t getFiltersb()   const = 0;
};

/** @class ObfFilterStatus
* @brief A container for the output of the OnboardFilter. Will contain a list of 
*        output objects for each of the filters run
*  
*  $Header$
*/

class ObfFilterStatus : public TObject  
{
public: 
    ObfFilterStatus();
    
    virtual ~ObfFilterStatus();

    // @brief Enumerate the keys for possible Filters we might have results for
    // @param GammaFilter - key for Gamma Filter results
    // @param CNOFilter   - key for CNO Filter results
    // @param MipFilter   - key for Mip Filter results
    enum FilterKeys {
        GammaFilter = 0,
        HipFilter   = 1,
        MipFilter   = 2,
        DgnFilter   = 3
    };
    
    /// clear lists, free pointers, etc., after read from / write to file
    void Clear(Option_t *option ="");
    void Fake( Int_t ievent, UInt_t rank, Float_t randNum ) ; // for tests
    Bool_t CompareInRange( const ObfFilterStatus &, const std::string & name = "" ) const ; // for tests
    void Print(Option_t *option="") const;
    /// Assignment operator for convenience
    ObfFilterStatus& operator=(const ObfFilterStatus& rhs);
    
    void initialize();

    // Add results method
//    void addFilterStatus(FilterKeys key, IObfStatus* status);
    void addFilterStatus(FilterKeys key, TObject* status);

    // Return results method
    const IObfStatus* getFilterStatus(FilterKeys key) const;

private:
    /// An array of points defining the trajectory
    TObjArray m_obfStatusCol;
    
    ClassDef(ObfFilterStatus,3) // Onboard Filter container class
};

class ObfGammaStatus : virtual public IObfStatus, public TObject
{
public:
    ObfGammaStatus() : m_status(0), m_id(0), m_sb(0), m_energy(0) {}
    ObfGammaStatus(UChar_t id, UInt_t status, UChar_t sb, UInt_t energy) : 
                   m_status(status), m_id(id), m_sb(sb), m_energy(energy) {}
    virtual ~ObfGammaStatus() {}

    // If msb is set below then event is to be vetoed
    UInt_t  getStatusWord() const;

    UInt_t  getVetoMask()   const;
    UInt_t  getVetoBit()    const;

    UChar_t getFilterId()   const {return m_id;}
    UChar_t getFiltersb()   const {return m_sb;}

    UInt_t  getEnergy()     const {return m_energy;}
    
    void Clear(Option_t *option ="");
    void Print(Option_t *option="") const;

private:
    UInt_t  m_status;
    UChar_t m_id;
    UChar_t m_sb;
    UInt_t  m_energy;

    ClassDef(ObfGammaStatus,4) // Gamma Filter output
};

class ObfHipStatus : virtual public IObfStatus, public TObject
{
public:
    ObfHipStatus() : m_status(0), m_id(0), m_sb(0) {}
    ObfHipStatus(UChar_t id, UInt_t status, UChar_t sb) : 
                   m_status(status), m_id(id), m_sb(sb) {}
    virtual ~ObfHipStatus() {}

    // If msb is set below then event is to be vetoed
    UInt_t  getStatusWord() const;

    UInt_t  getVetoMask()   const;
    UInt_t  getVetoBit()    const;

    UChar_t getFilterId()   const {return m_id;}
    UChar_t getFiltersb()   const {return m_sb;}
    
    void Clear(Option_t *option ="");
    void Print(Option_t *option="") const;

private:
    UInt_t  m_status;
    UChar_t m_id;
    UChar_t m_sb;

    ClassDef(ObfHipStatus,4) // HFC output
};

class ObfMipStatus : virtual public IObfStatus, public TObject
{
public:
    ObfMipStatus(): m_status(0), m_id(0), m_sb(0) {}
    ObfMipStatus(UChar_t id, UInt_t status, UChar_t sb) : 
                   m_status(status), m_id(id), m_sb(sb) {}
    virtual ~ObfMipStatus() {}

    // If msb is set below then event is to be vetoed
    UInt_t  getStatusWord() const;

    UInt_t  getVetoMask()   const;
    UInt_t  getVetoBit()    const;

    UChar_t getFilterId()   const {return m_id;}
    UChar_t getFiltersb()   const {return m_sb;}
    
    void Clear(Option_t *option ="");
    void Print(Option_t *option="") const;

private:
    UInt_t  m_status;
    UChar_t m_id;
    UChar_t m_sb;

    ClassDef(ObfMipStatus,4) // MIP filter output
};

class ObfDgnStatus : virtual public IObfStatus, public TObject
{
public:
    ObfDgnStatus() : m_status(0), m_id(0), m_sb(0) {}
    ObfDgnStatus(UChar_t id, UInt_t status, UChar_t sb)
     : m_status(status), m_id(id), m_sb(sb) {}
    virtual ~ObfDgnStatus() {}

    // If msb is set below then event is to be vetoed
    UInt_t  getStatusWord() const;

    UInt_t  getVetoMask()   const;
    UInt_t  getVetoBit()    const;

    UChar_t getFilterId()   const {return m_id;}
    UChar_t getFiltersb()   const {return m_sb;}
    
    void Clear(Option_t *option ="");
    void Print(Option_t *option="") const;

private:
    UInt_t  m_status;
    UChar_t m_id;
    UChar_t m_sb;

    ClassDef(ObfDgnStatus,4) // DFC output
};

#endif
