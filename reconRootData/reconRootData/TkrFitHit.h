#ifndef ROOT_TKRFITHIT_H
#define ROOT_TKRFITHIT_H

#include "TObject.h"
#include "TkrParams.h"
#include "TkrCovMat.h"
/** 
 * @class TkrFitHit
 *
 * @brief Root object containing fit track parameters and covariance matrix at one point
 *        This class should duplicate Event::TkrFitHit
 * 
 * Created 14-May-2002
 * 
 * @author Tracy Usher
 *
 */

class TkrFitHit : public TObject
{
public:
    //! This probably does not want to be 
    enum TYPE{MEAS, FIT, PRED, SMOOTH, UNKNOWN};

    TkrFitHit();

    TkrFitHit(TYPE type, const TkrParams& params, const TkrCovMat& covMat);

    TkrFitHit(const TkrFitHit& hit);

    virtual ~TkrFitHit() {}

    TkrFitHit& operator=(const TkrFitHit& old);

    //! Provide an initialize method (in the case of the default constructor)
    void     initialize(TYPE type, const TkrParams& params, const TkrCovMat& covMat);

    //! Group the methods for returning info here
    TYPE             getHitType()   const {return m_HitType;}
    const TkrParams& getTkrParams() const {return m_Params;}
    const TkrCovMat& getTkrCovMat() const {return m_CovMat;}

private:
    TYPE       m_HitType;     // Type of fit track parameters contained within
    TkrParams  m_Params;      // Track parameters at this point
    TkrCovMat  m_CovMat;      // Covariance matrix associated with these parameters

    ClassDef(TkrFitHit,1)
};

#endif
