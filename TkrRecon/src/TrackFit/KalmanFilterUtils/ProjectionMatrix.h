/**
 * @class ProjectionMatrix
 *
 * @brief Definition of a Kalman Filter transport matrix
 *
 * @author Tracy Usher
 *
 * $Header$
 */

#ifndef ProjectionMatrix_h
#define ProjectionMatrix_h

#include "KalmanFilterDefs.h"

class ProjectionMatrix 
{
public:

    virtual KFmatrix H(int i) = 0;
    virtual KFmatrix operator()(const int &i) = 0;

};


#endif
