//  version 0.1
//
//  Energy of a calorimeter cluster
//
#ifndef CALLOGENE_H
#define CALLOGENE_H

// Root Includes
#include "TObject.h"
#include "TVector3.h"
// Local Cal forward declaration
#include "digiRootData/LogID.h"

// =========================================================
//      CalLogEne
// =========================================================

//! Root object that contains the response of a single CsI log
/*! Contains:
    1.) LogID
    2.) Response from the Minus end
    2.) Response from the Plus end */
class CalLogEne : public TObject
{

private:

    //! response from either end of the log
    Float_t  m_eMinus;    
    Float_t  m_ePlus; 
    
    //! Pointer to the identifier for this log
    LogID *m_id;

  
public:

    //! default constructor
    CalLogEne();
    //! constructor
    CalLogEne(UInt_t layer, UInt_t column, LogID::CALAxes xy);

    //! destructor
    virtual ~CalLogEne();

    void Clean();

    //! access the LogID 
    LogID* getId () { return m_id; };
    //! set the LogID
    inline void setId(LogID *id) { m_id = id; };

    //!  access response on the MINUS end of the log
    Float_t getEminus() { return m_eMinus; };
    //! set the response on the MINUS end of the log
    void seteMinus(Float_t eM) { m_eMinus = eM; };

    //! access the response on the PLUS end of the log
    Float_t getEplus() { return m_ePlus; };
    //! set the response on the PLUS end of the log
    void setePlus(Float_t eP) { m_ePlus = eP; };

    ClassDef(CalLogEne, 1)
};

#endif
