/**
 * @class GlastTransportMatrix
 *
 * @brief Definition of a Kalman Filter transport matrix
 *
 * @author Tracy Usher
 *
 * $Header$
 */

#ifndef GlastTransportMatrix_h
#define GlastTransportMatrix_h

#include "src/TrackFit/KalmanFilterUtils/TransportMatrix.h"
#include <vector>

class GlastTransportMatrix : public TransportMatrix
{
public:

    // Constructor 
    GlastTransportMatrix(std::vector<double> zCoords);
   ~GlastTransportMatrix() {};

    virtual KFmatrix F(int i, int j);
    virtual KFmatrix operator()(const int &i, const int &j);

private:
    std::vector<double> m_zCoords;
};


#endif