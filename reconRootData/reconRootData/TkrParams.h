#ifndef ROOT_TKRPARAMS_H
#define ROOT_TKRPARAMS_H

#include "TObject.h"
//#include "TVectorD.h"
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

class TkrParams : public TObject
{
public:

    TkrParams();

    TkrParams(Double_t ax, Double_t sx, Double_t ay, Double_t sy);

    TkrParams(const TkrParams& params);

    virtual ~TkrParams() {};

    /// resets all member variables to default values
    void Clear(Option_t *option="");
	
    void Print(Option_t *option="") const;

    TkrParams& operator=(const TkrParams& old);

    /// Provide an initialization method (in case of default constructor)
    void initialize(Double_t ax, Double_t sx, Double_t ay, Double_t sy);
    void initialize(const TkrParams& params);

    const Double_t& TkrParams::operator() (const int &element) const;
    Double_t& TkrParams::operator() (const int &element);

    /// NOTE index starts with [0]. In case of index error returns 0
    Double_t& TkrParams::operator[] (const int &element);
    const Double_t& TkrParams::operator[] (const int &element) const;

    /// operator overload for TkrParams+TkrParams
    const TkrParams operator +(const TkrParams& A) const;
	      
    /// operator overload for TkrParams-TkrParams
    const TkrParams operator -(const TkrParams& A) const;
		
    /// operator overload for TkrParams*TkrParams
    const Double_t operator *(const TkrParams& A) const;

    inline Double_t getXPos()   const { return m_xpos; }
    inline Double_t getYPos()   const { return m_ypos; }
    inline Double_t getXSlope() const { return m_xslope; }
    inline Double_t getYSlope() const { return m_yslope; }

private:

    Double_t m_xpos;
    Double_t m_xslope;
    Double_t m_ypos;
    Double_t m_yslope;

    ClassDef(TkrParams,2)

};
#endif
