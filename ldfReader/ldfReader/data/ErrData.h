#ifndef ldfReader_ErrData_H
#define ldfReader_ErrData_H

#include <vector>

#include "EventSummaryCommon.h"

namespace ldfReader {

    /** @class ErrData
      * @brief Local storage of error data
      * $Header$
    */
    class ErrData {
    public:

        ErrData() { clear(); };
        ErrData(const ErrData& error) { m_summary = error.m_summary; };
        ErrData(const EventSummaryCommon &summary) { m_summary = summary; }
        ~ErrData() { clear(); };

        void clear() { m_summary.setSummary(0); };

        const EventSummaryCommon& summary() const { return m_summary; };
        void setSummary(unsigned summary) { m_summary.setSummary(summary);};

    private:

        // Store the event sequence number for this contribution
        EventSummaryCommon m_summary;

    };
} // end namespace
#endif
