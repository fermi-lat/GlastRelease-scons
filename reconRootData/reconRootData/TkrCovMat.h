#ifndef ROOT_TKRCOVMAT_H
#define ROOT_TKRCOVMAT_H

#include "TMatrixD.h"
/** 
 * @class TkrCovMat
 * @brief Root object containing track covariance matrix at a point
 *
 * This class should duplicate Event::TkrFitMatrix TDS class
 * 
 * @author Tracy Usher
 *
 * $Header$
 */

class TkrCovMat : public TMatrixD
{
public:

    TkrCovMat() : TMatrixD(4,4) { Zero(); };

    TkrCovMat(const TkrCovMat& covMat);

    TkrCovMat( Double_t e_00, Double_t e_01, Double_t e_02, Double_t e_03,
        Double_t e_10, Double_t e_11, Double_t e_12, Double_t e_13,
        Double_t e_20, Double_t e_21, Double_t e_22, Double_t e_23,
        Double_t e_30, Double_t e_31, Double_t e_32, Double_t e_33);


    virtual ~TkrCovMat() {}

    void Print(Option_t *option="") const;

    TkrCovMat& operator=(const TkrCovMat& covMat);

    void     initialize(const TkrCovMat& covMat);

    void initialize(Double_t e_00, Double_t e_01, Double_t e_02, Double_t e_03,
        Double_t e_10, Double_t e_11, Double_t e_12, Double_t e_13,
        Double_t e_20, Double_t e_21, Double_t e_22, Double_t e_23,
        Double_t e_30, Double_t e_31, Double_t e_32, Double_t e_33);						 


    /// Methods to return errors 
    inline Double_t getCovX0X0() const {return GetNoElements() == 16 ? operator()(0,0) : -9999.;}
    inline Double_t getCovX0Sx() const {return GetNoElements() == 16 ? operator()(0,1) : -9999.;}
    inline Double_t getCovX0Y0() const {return GetNoElements() == 16 ? operator()(0,2) : -9999.;}
    inline Double_t getCovX0Sy() const {return GetNoElements() == 16 ? operator()(0,3) : -9999.;}

    inline Double_t getCovSxX0() const {return GetNoElements() == 16 ? operator()(1,0) : -9999.;}
    inline Double_t getCovSxSx() const {return GetNoElements() == 16 ? operator()(1,1) : -9999.;}
    inline Double_t getCovSxY0() const {return GetNoElements() == 16 ? operator()(1,2) : -9999.;}
    inline Double_t getCovSxSy() const {return GetNoElements() == 16 ? operator()(1,3) : -9999.;}

    inline Double_t getCovY0X0() const {return GetNoElements() == 16 ? operator()(2,0) : -9999.;}
    inline Double_t getCovY0Sx() const {return GetNoElements() == 16 ? operator()(2,1) : -9999.;}
    inline Double_t getCovY0Y0() const {return GetNoElements() == 16 ? operator()(2,2) : -9999.;}
    inline Double_t getCovY0Sy() const {return GetNoElements() == 16 ? operator()(2,3) : -9999.;}

    inline Double_t getCovSyX0() const {return GetNoElements() == 16 ? operator()(3,0) : -9999.;}
    inline Double_t getCovSySx() const {return GetNoElements() == 16 ? operator()(3,1) : -9999.;}
    inline Double_t getCovSyY0() const {return GetNoElements() == 16 ? operator()(3,2) : -9999.;}
    inline Double_t getCovSySy() const {return GetNoElements() == 16 ? operator()(3,3) : -9999.;}

private:

    ClassDef(TkrCovMat,1)
};
#endif
