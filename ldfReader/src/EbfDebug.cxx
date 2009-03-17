// $Header$
#include "EbfDebug.h"

int EbfDebug::setDebug(int debugOn) {
    bool old = m_debug;
    m_debug = debugOn;
    return old;
}

int EbfDebug::m_debug;
