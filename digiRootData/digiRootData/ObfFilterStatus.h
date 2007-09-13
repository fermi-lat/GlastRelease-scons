#ifndef obfFilterStatus_H
#define obfFilterStatus_H

#include "TObject.h"
#include "TObjArray.h"

/** @class IObfStatus
* @brief Interface/base class which describes access to the results of the various 
*        onboard filters. Currently, all filters output a 32 bit status word but then
*        can also include some number of other items. This interface class deals 
*        entirely with the status bits
*
*/
class IObfStatus : public TObject
{
public:
    virtual ~IObfStatus() {}

    /// TObject functions
    virtual void Clear(Option_t *option ="") = 0;
    virtual void Print(Option_t *option="") const = 0;

    /// OBF functions
    virtual UInt_t getStatusHi() const = 0;
    virtual UInt_t getStatusLo() const = 0;
    virtual UInt_t getStatus32() const = 0;
    virtual UInt_t getVetoMask() const = 0;
    virtual UInt_t getVetoBit()  const = 0;
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
        HFCFilter   = 1,
        MipFilter   = 2,
        DFCFilter   = 3
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
    void addFilterStatus(FilterKeys key, IObfStatus* status);

    // Return results method
    const IObfStatus* getFilterStatus(FilterKeys key);

private:
    /// An array of points defining the trajectory
    TObjArray m_obfStatusCol;
    
    ClassDef(ObfFilterStatus,1) // Onboard Filter container class
};

class ObfGammaStatus : virtual public IObfStatus
{
public:
    ObfGammaStatus()              : m_status(0)      {}
    ObfGammaStatus(UInt_t status) : m_status(status) {}
    virtual ~ObfGammaStatus() {}

    // This for backward compatibility...
    UInt_t getStatusHi() const; 
    UInt_t getStatusLo() const;
    // If msb is set below then event is to be vetoed
    UInt_t getStatus32() const;

    UInt_t getVetoMask() const;
    UInt_t getVetoBit()  const;
    
    void Clear(Option_t *option ="");
    void Print(Option_t *option="") const;

private:
    UInt_t m_status;

    ClassDef(ObfGammaStatus,1) // Gamma Filter output
};

class ObfHFCStatus : virtual public IObfStatus
{
public:
    ObfHFCStatus()              : m_status(0)      {}
    ObfHFCStatus(UInt_t status) : m_status(status) {}
    virtual ~ObfHFCStatus() {}

    // This for backward compatibility...
    UInt_t getStatusHi() const; 
    UInt_t getStatusLo() const;
    // If msb is set below then event is to be vetoed
    UInt_t getStatus32() const;

    UInt_t getVetoMask() const;
    UInt_t getVetoBit()  const;
    
    void Print(Option_t *option="") const;

private:
    UInt_t m_status;

    ClassDef(ObfHFCStatus,1) // HFC output
};

class ObfMipStatus : virtual public IObfStatus
{
public:
    ObfMipStatus()              : m_status(0)      {}
    ObfMipStatus(UInt_t status) : m_status(status) {}
    virtual ~ObfMipStatus() {}

    // This for backward compatibility...
    UInt_t getStatusHi() const; 
    UInt_t getStatusLo() const;
    // If msb is set below then event is to be vetoed
    UInt_t getStatus32() const;

    UInt_t getVetoMask() const;
    UInt_t getVetoBit()  const;
    
    void Print(Option_t *option="") const;

private:
    UInt_t m_status;

    ClassDef(ObfMipStatus,1) // MIP filter output
};

class ObfDFCStatus : virtual public IObfStatus
{
public:
    ObfDFCStatus()              : m_status(0)      {}
    ObfDFCStatus(UInt_t status) : m_status(status) {}
    virtual ~ObfDFCStatus() {}

    // This for backward compatibility...
    UInt_t getStatusHi() const; 
    UInt_t getStatusLo() const;
    // If msb is set below then event is to be vetoed
    UInt_t getStatus32() const;

    UInt_t getVetoMask() const;
    UInt_t getVetoBit()  const;
    
    void Print(Option_t *option="") const;

private:
    UInt_t m_status;

    ClassDef(ObfDFCStatus,1) // DFC output
};

#endif
