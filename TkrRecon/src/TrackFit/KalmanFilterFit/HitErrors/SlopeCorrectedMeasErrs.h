/**
 * @class SlopeCorrectedMeasErrs
 *
 * @brief This class assigns errors to the measured coordinates in a given plane
 *        SlopeCorrectedMeasErrs implements a version which attempts to provide a detailed measurement 
 *        error based on cluster width, incoming slope, etc. 
 *        ** Experimental version **
 *
 * @author Tracy Usher (editor) from version implemented by Leon Rochester (due to Bill Atwood)
 *
 * $Header$
 */

#ifndef SlopeCorrectedMeasErrs_h
#define SlopeCorrectedMeasErrs_h

#include "IComputeMeasErrors.h"
#include "src/Track/TkrControl.h"
#include "TkrUtil/ITkrGeometrySvc.h"

class SlopeCorrectedMeasErrs : public IComputeMeasErrors
{
public:

    // Constructor needs the matrices that transform state vector, covariance matrix
    SlopeCorrectedMeasErrs(ITkrGeometrySvc* tkrGeo);
    virtual ~SlopeCorrectedMeasErrs() {};


    TkrCovMatrix computeMeasErrs(const Event::TkrTrackParams& newPars, 
                                 const TkrCovMatrix&          oldCovMat, 
                                 const Event::TkrCluster&     cluster);

private:
    double getError(double strips, double slope) const;

    ITkrGeometrySvc* m_tkrGeo;
    TkrControl*      m_control;
};


#endif
