#include "reconRootData/TkrParams.h"

ClassImp(TkrParams)

void TkrParams::initialize(Double_t ax, Double_t sx, Double_t ay, Double_t sy)
{
    if (GetNoElements() != 4) ResizeTo(4);

    Zero();

    operator()(0) = ax;
    operator()(1) = sx;
    operator()(2) = ay;
    operator()(3) = sy;
}  
