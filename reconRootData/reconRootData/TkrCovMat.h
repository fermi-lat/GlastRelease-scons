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

    TkrCovMat(Double_t sigXpos, Double_t sigXslp, Double_t sigXps, 
              Double_t sigYpos, Double_t sigYslp, Double_t sigYps); 

    virtual ~TkrCovMat() {}

    void Print(Option_t *option="") const;

    TkrCovMat& operator=(const TkrCovMat& covMat);

    //! Provide an initialization method (in case of default constructor)
    void     initialize(Double_t sigXpos, Double_t sigXslp, Double_t sigXps, 
                        Double_t sigYpos, Double_t sigYslp, Double_t sigYps);
    void     initialize(const TkrCovMat& covMat);

    //! Methods to return errors 
    //! (Are these redundant with the inheritance of TMatrixD?)
    Double_t getCovX0X0() const {return GetNoElements() == 16 ? operator()(0,0) : -9999.;}
    Double_t getCovSxSx() const {return GetNoElements() == 16 ? operator()(1,1) : -9999.;}
    Double_t getCovX0Sx() const {return GetNoElements() == 16 ? operator()(0,1) : -9999.;}
    Double_t getCovY0Y0() const {return GetNoElements() == 16 ? operator()(2,2) : -9999.;}
    Double_t getCovSySy() const {return GetNoElements() == 16 ? operator()(3,3) : -9999.;}
    Double_t getCovY0Sy() const {return GetNoElements() == 16 ? operator()(2,3) : -9999.;}

private:

    ClassDef(TkrCovMat,1)
};
#endif
