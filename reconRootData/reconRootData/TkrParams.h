#ifndef ROOT_TKRPARAMS_H
#define ROOT_TKRPARAMS_H

#include "TObject.h"
#include "TVectorD.h"
/** 
 * @class TkrParams
 * @brief Root object containing fit track parameters for a given plane
 * 
 * This class should be a duplicate of Event::TkrFitPar
 * 
 * @author Tracy Usher
 *
 * $Header$
 */

class TkrParams : public TVectorD
{
public:

    TkrParams() {};

    TkrParams(Double_t ax, Double_t sx, Double_t ay, Double_t sy) : TVectorD(0,3,ax,sx,ay,sy,"END") 
    {}

    TkrParams(const TkrParams& params);

    virtual ~TkrParams() {};

    TkrParams& operator=(const TkrParams& old);

    //! Provide an initialization method (in case of default constructor)
    void     initialize(Double_t ax, Double_t sx, Double_t ay, Double_t sy);
    void     initialize(const TkrParams& params);

    //! Get methods
    Double_t getXPos()   const {return GetNoElements() == 4 ? operator()(0) : -9999.;}
    Double_t getYPos()   const {return GetNoElements() == 4 ? operator()(2) : -9999.;}
    Double_t getXSlope() const {return GetNoElements() == 4 ? operator()(1) : -9999.;}
    Double_t getYSlope() const {return GetNoElements() == 4 ? operator()(3) : -9999.;}

    ClassDef(TkrParams,1)
private:
};
#endif
