/**
 * @class TransportMatrix
 *
 * @brief Definition of a Kalman Filter transport matrix
 *
 * @author Tracy Usher
 *
 * $Header$
 */

#ifndef TransportMatrix_h
#define TransportMatrix_h

#include "KalmanFilterDefs.h"

class TransportMatrix 
{
public:

    virtual KFmatrix F(int i, int j) = 0;
    virtual KFmatrix operator()(const int &i, const int &j) = 0;
};


#endif