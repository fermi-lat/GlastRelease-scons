#ifndef ldfReader_OswData_H
#define ldfReader_OswData_H

#include <vector>

#include "EventSummaryCommon.h"

namespace ldfReader {

    /** @class OswData
      * @brief Local storage of OSW data
      * $Header$
    */
    class OswData {
    public:

        OswData() { clear(); };
        OswData(const OswData& osw) { 
            m_summary = osw.m_summary; 
            m_lenInBytes = osw.m_lenInBytes; };
        OswData(const EventSummaryCommon& summary) { m_summary = summary; };
        ~OswData() { clear(); };

        void clear() { 
            m_summary.setSummary(0); 
            m_lenInBytes = 0; };

        const EventSummaryCommon& summary() const { return m_summary; };
        void initSummary(unsigned summary) { m_summary.setSummary(summary);};

        void initLength(unsigned long length) { m_lenInBytes = length; };
        unsigned long lenInBytes() const { return m_lenInBytes; };

    private:

        // Store the event sequence number for this contribution
        EventSummaryCommon m_summary;

        unsigned long m_lenInBytes;

    };
} // end namespace
#endif
