#ifndef TKRFITHIT_H
#define TKRFITHIT_H

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

    //! constructors
    TkrFitHit();

    TkrFitHit(TYPE type, const TkrParams& params, const TkrCovMat& covMat) :
              m_HitType(type),
              m_Params(params),
              m_CovMat(covMat)
              {}

    TkrFitHit(const TkrFitHit& hit);

    //! destructor
    virtual ~TkrFitHit() {}

    //! Provide an initialize method (in the case of the default constructor)
    void     initialize(TYPE type, const TkrParams& params, const TkrCovMat& covMat);

    //! Group the methods for returning info here
    TYPE             getHitType()   const {return m_HitType;}
    const TkrParams& getTkrParams() const {return m_Params;}
    const TkrCovMat& getTkrCovMat() const {return m_CovMat;}

    ClassDef(TkrFitHit,1)
private:
    TYPE       m_HitType;     // Type of fit track parameters contained within
    TkrParams  m_Params;      // Track parameters at this point
    TkrCovMat  m_CovMat;      // Covariance matrix associated with these parameters
};

#endif
