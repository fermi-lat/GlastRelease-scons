#ifndef ROOT_Tem_H
#define ROOT_Tem_H

#include "TObject.h"
#include "ErrorData.h"


/** @class Tem
 * @brief This stores TEM data.
 *
 * The Tem branch contains:
 *
 * $Header$
 */

class Tem: public TObject {


public:

    Tem();
    Tem(UShort_t t, const ErrorData &err);
    Tem(const Tem& t); 
    virtual ~Tem() { };
    void init(UShort_t t, const ErrorData &err);

    void Clear(Option_t *option="");
 
    void Print(Option_t *option="") const;



    /** @ingroup TemGroup */
    /*@{*/
    /// ID of the TEM, which does correspond to Gleam's tower numbering
    UShort_t getTowerId() const { return m_tower; };
    /// Retrieve the error contribution for this TEM, which may be empty
    const ErrorData& getError() const { return m_error; };
    /*@}*/

private:

    UShort_t m_tower;
    ErrorData m_error;

    ClassDef(Tem,1) // Storage for TEM 
}; 
 
#endif


