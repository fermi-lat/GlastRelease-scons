/**
 * @class IKalmanFilterMatrix
 *
 * @brief Implementation Process Noise matrix for the Kalman Filter
 *
 * @author Tracy Usher
 *
 * $Header$
 */

#ifndef IKalmanFilterMatrix_h
#define IKalmanFilterMatrix_h

#include "KalmanFilterDefs.h"

class KalmanFilterInit;

class IKalmanFilterMatrix 
{
public:
    virtual void     accept(const KalmanFilterInit& init) = 0;

    virtual KFmatrix operator()(const int &i) = 0;
    virtual KFmatrix operator()(const int &i, const int &j) = 0;
    virtual KFmatrix operator()(const KFvector& stateVec, const int &i, const int &j) = 0;
};


#endif
