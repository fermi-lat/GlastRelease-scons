/**
 * @class StdTransportMatrix
 *
 * @brief Definition of a Kalman Filter transport matrix
 *
 * @author Tracy Usher
 *
 * $Header$
 */

#include "StdTransportMatrix.h"

KFmatrix& StdTransportMatrix::operator ()(const double& deltaZ)
{
    m_F(1,2) = deltaZ;
    m_F(3,4) = deltaZ;

    return m_F;
}

