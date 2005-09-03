/** @file TkrTruncatedPlane.h
* @author Leon Rochester
*
* $Header$

*/

#ifndef TkrTruncatedPlane_H
#define TkrTruncatedPlane_H

#include <vector>

static const CLID& CLID_TkrTruncatedPlane = InterfaceID("TkrTruncatedPlane",  0, 0);


typedef std::vector<int>   intVector;
typedef std::vector<float> floatVector;

namespace Event {
/**
* @class TkrTruncatedPlane
*
* @brief Keeps track of hit truncation information for a single plane
*
* TkrTruncatedPlane 
*
*/

    class TkrTruncatedPlane {

    public:
        TkrTruncatedPlane(int status=0, 
            intVector   stripCount=intVector(0), 
            intVector   stripNumber=intVector(0),
            floatVector localX=floatVector(0),
            float planeZ=0) : m_status(status), m_stripCount(stripCount), 
            m_stripNumber(stripNumber), m_localX(localX), m_planeZ(planeZ) {}

            ~TkrTruncatedPlane() {}

            enum truncationType
            {
                RC       = 0x01,
                CC       = 0x02,
                RCOVER   = 0x04,
                CCOVER   = 0x08,
                shiftEND = 4,
                maskEND  = 0xF,
                END0SET  = maskEND,
                END1SET  = (END0SET<<shiftEND),
                RC0SET   = (RC | RCOVER),
                CC0SET   = (CC | CCOVER),
                RC1SET   = (RC0SET<<shiftEND),
                CC1SET   = (CC0SET<<shiftEND),
                RCSET    = (RC0SET | RC1SET),
                CCSET    = (CC0SET | CC1SET)
            };

            /// adds a bit or bit pattern to the status word
            void setStatusBit(int end, TkrTruncatedPlane::truncationType type) {
                m_status = (m_status |(end==0 ? type : type<<shiftEND));
            }
            /// sets the whole status word
            void setStatus(int status) { m_status = status; }
            /// retrieves the status word
            const int          getStatus() const      { return m_status; }
            /// retrieves a reference to the numStrips vector
            const intVector&   getStripCount() const   { return m_stripCount; }
            /// retrieves a reference to the stripNumber vector
            const intVector&   getStripNumber() const { return m_stripNumber; }
            /// retrieves a reference to the stripPosition vector
            const floatVector& getLocalX() const      { return m_localX; }
            /// retrieves the global z of the plane
            const float        getPlaneZ() const      { return m_planeZ; }
            /// true if there are truncation planes in the event
            const bool         isTruncated() const    { return  (m_status != 0); }
            /// sets the status bit in a non-member status word
            static void addStatusBit(int& status, int end, TkrTruncatedPlane::truncationType type) {
                status = (status | (end==0 ? type : type<<shiftEND));
            }

    private:
        /// status of plane, see bit definitions
        int m_status;
        /// number of strips on each end
        intVector m_stripCount;
        /// strip numbers: highest low, lowest high, highest high and splitpoint
        intVector m_stripNumber;
        /// local coordinates of the strip numbers above
        floatVector m_localX;
        /// z coordinate of the plane
        float m_planeZ;
    };
}

#endif
