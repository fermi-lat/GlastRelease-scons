/**
 * @class StdProjectionMatrix
 *
 * @brief Definition of a Kalman Filter projection matrix. This is used to projected out the measured
 *        coordinates (and errors) from the state vector (covariance matrix) during the track fits
 *
 * @author Tracy Usher
 *
 * $Header$
 */

#include "StdProjectionMatrix.h"
#include "src/TrackFit/KalmanFilterFit/KalmanFilterInit.h"
#include "idents/TkrId.h"


StdProjectionMatrix::StdProjectionMatrix() : m_none(1,4), m_projX(1,4), m_projY(1,4)
{
    m_projX(1,1) = 1;
    m_projY(1,3) = 1;

    return;
}

KFmatrix& StdProjectionMatrix::operator()(const idents::TkrId &id)
{
    if (id.getView() == idents::TkrId::eMeasureX) return m_projX;

    return m_projY;
}
