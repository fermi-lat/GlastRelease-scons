#include "reconRootData/TkrParams.h"

ClassImp(TkrParams)

TkrParams::TkrParams(const TkrParams& params)
{
    initialize(params);
}

TkrParams& TkrParams::operator=(const TkrParams& params)
{
    initialize(params);

    return *this;
}

void TkrParams::initialize(Double_t ax, Double_t sx, Double_t ay, Double_t sy)
{
    if (GetNoElements() != 4) ResizeTo(4);

    operator()(0) = ax;
    operator()(1) = sx;
    operator()(2) = ay;
    operator()(3) = sy;
}  

void TkrParams::initialize(const TkrParams& params)
{
    ResizeTo(params);

    operator()(0) = params.getXPos();
    operator()(1) = params.getXSlope();
    operator()(2) = params.getYPos();
    operator()(3) = params.getYSlope();
}  
