/**
 * @class ProjectionMatrix
 *
 * @brief Definition of a Kalman Filter projection matrix
 *
 * @author Tracy Usher
 *
 * $Header$
 */

#ifndef GlastProjectionMatrix_h
#define GlastProjectionMatrix_h

#include "src/TrackFit/KalmanFilterUtils/ProjectionMatrix.h"
#include "CLHEP/Matrix/Matrix.h"
#include <vector>

class GlastProjectionMatrix : public ProjectionMatrix
{
public:

    // Constructor 
    GlastProjectionMatrix(std::vector<int> projection);
    ~GlastProjectionMatrix() {};

    KFmatrix H(int i);
    KFmatrix operator()(const int &i);

private:
    std::vector<int> m_projection;
};


#endif