#include "reconRootData/TkrCovMat.h"

ClassImp(TkrCovMat)

TkrCovMat::TkrCovMat(Double_t sigXpos, Double_t sigXslp, Double_t sigXps, 
                     Double_t sigYpos, Double_t sigYslp, Double_t sigYps) :
           TMatrixD(4,4)
{
    initialize(sigXpos, sigXslp, sigXps, sigYpos, sigYslp, sigYps);
}



void TkrCovMat::initialize(Double_t sigXpos, Double_t sigXslp, Double_t sigXps, 
                           Double_t sigYpos, Double_t sigYslp, Double_t sigYps)
{
    if (GetNrows() != 4 || GetNcols() != 4) ResizeTo(4,4);

    Zero();
    
    operator()(0,0) = sigXpos;
    operator()(1,1) = sigXslp;
    operator()(0,1) = sigXps;
    operator()(1,0) = sigXps;

    operator()(2,2) = sigYpos;
    operator()(3,3) = sigYslp;
    operator()(2,3) = sigYps;
    operator()(3,2) = sigYps;
}
