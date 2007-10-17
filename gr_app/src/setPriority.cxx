#ifdef WIN32
#include "windows.h"
#include "winbase.h"
#endif
    void setPriority()
    {
#ifdef WIN32
        ::SetThreadPriority(::GetCurrentThread(),THREAD_PRIORITY_BELOW_NORMAL);
#endif
}