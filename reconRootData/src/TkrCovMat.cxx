#include "reconRootData/TkrCovMat.h"
#include <iostream>

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

TkrCovMat::TkrCovMat(double e_00, double e_01, double e_02, double e_03,
                     double e_10, double e_11, double e_12, double e_13,
                     double e_20, double e_21, double e_22, double e_23,
                     double e_30, double e_31, double e_32, double e_33) 
{
    initialize( e_00, e_01, e_02, e_03, e_10, e_11, e_12, e_13,
        e_20, e_21, e_22, e_23, e_30, e_31, e_32, e_33);
}


void TkrCovMat::initialize(double e_00, double e_01, double e_02, double e_03,
                           double e_10, double e_11, double e_12, double e_13,
                           double e_20, double e_21, double e_22, double e_23,
                           double e_30, double e_31, double e_32, double e_33) 
{
    if (GetNrows() != 4 || GetNcols() != 4) ResizeTo(4,4);

    Zero();
    
    operator()(0,0) = e_00;
    operator()(0,1) = e_01;
    operator()(0,2) = e_02;
    operator()(0,3) = e_03;
    
    operator()(1,0) = e_10;
    operator()(1,1) = e_11;
    operator()(1,2) = e_12;
    operator()(1,3) = e_13;

    operator()(2,0) = e_20;
    operator()(2,1) = e_21;
    operator()(2,2) = e_22;
    operator()(2,3) = e_23;

    operator()(3,0) = e_30;
    operator()(3,1) = e_31;
    operator()(3,2) = e_32;
    operator()(3,3) = e_33;

}


void TkrCovMat::initialize(const TkrCovMat& covMat)
{
    ResizeTo(covMat);
    
    operator()(0,0) = covMat.getCovX0X0();
    operator()(0,1) = covMat.getCovX0Sx();
    operator()(0,2) = covMat.getCovX0Y0();
    operator()(0,3) = covMat.getCovX0Sy();

    operator()(1,0) = covMat.getCovSxX0();
    operator()(1,1) = covMat.getCovSxSx();
    operator()(1,2) = covMat.getCovSxY0();
    operator()(1,3) = covMat.getCovSxSy();

    operator()(2,0) = covMat.getCovY0X0();
    operator()(2,1) = covMat.getCovY0Sx();
    operator()(2,2) = covMat.getCovY0Y0();
    operator()(2,3) = covMat.getCovY0Sy();

    operator()(3,0) = covMat.getCovSyX0();
    operator()(3,1) = covMat.getCovSySx();
    operator()(3,2) = covMat.getCovSyY0();
    operator()(3,3) = covMat.getCovSySy();
}



void TkrCovMat::Print(Option_t *option) const {
    TObject::Print(option);
    using namespace std;
    cout << operator()(0,0) << " " <<  operator()(0,1) << " "
        << operator()(0,2) << " " << operator()(0,3) << std::endl;
    cout << operator()(1,0) << " " << operator()(1,1) << " "
        << operator()(1,2) << " " << operator()(1,3) << std::endl;
    cout << operator()(2,0) << " " << operator()(2,1) << " "
        << operator()(2,2) << " " << operator()(2,3) << std::endl;
    cout << operator()(3,0) << " " << operator()(3,1) << " "
        << operator()(3,2) << " " << operator()(3,3) << std::endl;

}