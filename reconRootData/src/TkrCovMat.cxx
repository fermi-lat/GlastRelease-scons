#include "reconRootData/TkrCovMat.h"

ClassImp(TkrCovMat)

TkrCovMat::TkrCovMat(const TkrCovMat& covMat)
{
    initialize(covMat);
}

TkrCovMat& TkrCovMat::operator=(const TkrCovMat& covMat)
{
    initialize(covMat);

    return *this;
}

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

void TkrCovMat::initialize(const TkrCovMat& covMat)
{
    ResizeTo(covMat);
    
    operator()(0,0) = covMat.getCovX0X0();
    operator()(1,1) = covMat.getCovSxSx();
    operator()(0,1) = covMat.getCovX0Sx();
    operator()(1,0) = covMat.getCovX0Sx();
    
    operator()(0,0) = covMat.getCovY0Y0();
    operator()(1,1) = covMat.getCovSySy();
    operator()(0,1) = covMat.getCovY0Sy();
    operator()(1,0) = covMat.getCovY0Sy();
}
