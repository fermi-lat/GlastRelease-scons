#ifndef ROOT_ErrorData_H
#define ROOT_ErrorData_H

#include "TObject.h"
#include "enums/DetectorConstants.h"


/** @class ErrorData
 * @brief This stores Error data for one TEM.
 *
 * Read more about the error contribution at:
 * http://www-glast.slac.stanford.edu/IntegrationTest/ONLINE/docs/TEM.pdf 
 *
 * $Header$
 */

class ErrorData: public TObject {


public:

    ErrorData();
    ErrorData(UShort_t cal, UShort_t tkr, Bool_t phs, Bool_t tmo);
    ErrorData(const ErrorData& err); 
    virtual ~ErrorData() { };

    void init(UShort_t cal, UShort_t tkr, Bool_t phs, Bool_t tmo);
    void setTkrFifoFull(UChar_t* tkrFifoFullCol);
    void setTkrFifoFull(UInt_t gtcc, UChar_t val=1);

    void Clear(Option_t *option="");
 
    void Print(Option_t *option="") const;


    /** @defgroup ErrorDataGroup ErrorData End-User Interface */
    /** @ingroup ErrorDataGroup */
    /*@{*/
        UShort_t getCal() const { return m_cal; }
        UShort_t getTkr() const { return m_tkr; }
        Bool_t getPhs() const { return m_phs; }
        Bool_t getTmo() const { return m_tmo; } 
        const UChar_t* getTkrFifoFullCol() const { return m_tkrFifoFullCol; }
        Int_t getTkrFifoFull(UInt_t gtcc) const { 
            if (gtcc < enums::numGtcc) 
                return ((Int_t)m_tkrFifoFullCol[gtcc]);
            else
                return -1;
        }
    /*@}*/

    private:

        UShort_t m_cal;
        UShort_t m_tkr;
        Bool_t m_phs;
        Bool_t m_tmo;
        UChar_t m_tkrFifoFullCol[8];

    ClassDef(ErrorData,2) // Storage for TEM 
}; 
 
#endif


