#ifndef ROOT_ErrorData_H
#define ROOT_ErrorData_H

#include "TObject.h"


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

    void Clear(Option_t *option="");
 
    void Print(Option_t *option="") const;


    /** @defgroup ErrorDataGroup ErrorData End-User Interface */
    /** @ingroup ErrorDataGroup */
    /*@{*/
        UShort_t getCal() const { return m_cal; }
        UShort_t getTkr() const { return m_tkr; }
        Bool_t getPhs() const { return m_phs; }
        Bool_t getTmo() const { return m_tmo; } 
    /*@}*/

    private:

        UShort_t m_cal;
        UShort_t m_tkr;
        Bool_t m_phs;
        Bool_t m_tmo;

    ClassDef(ErrorData,1) // Storage for TEM 
}; 
 
#endif


