// $Header$
#include "EbfDebug.h"

bool EbfDebug::setDebug(bool debugOn) {
    bool old = m_debug;
    m_debug = debugOn;
    return old;
}
// Do this one inline
// bool EbfDebuf::getDebug() {return m_debug;}

bool EbfDebug::m_debug;
