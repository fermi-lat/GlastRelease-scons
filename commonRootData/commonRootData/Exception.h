#ifndef ROOT_EXCEPTION_H
#define ROOT_EXCEPTION_H

#include <stdio.h>

class BaseException {
protected:
    void PrintErrMsg() { 
        printf(m_szErrMsg);
        printf("\n");
    };
public:
    char *m_szErrMsg;
    
    BaseException() { m_szErrMsg = "no information available"; };
    
    void debugPrint() const { 
        printf("EXCEPTION: \n");
    };
};


class IndexOutOfBoundsException : public BaseException {
protected:
    void PrintErrMsg() { printf(m_szErrMsg, m_nIndex, m_nLowerBound, m_nUpperBound); };
    
public:
    IndexOutOfBoundsException(long nIndex, long nLower, long nUpper) { 
        m_nIndex = nIndex;
        m_nLowerBound = nLower; 
        m_nUpperBound = nUpper;
        m_szErrMsg = "Index, %d, out of bounds [%d, %d]";
    };
    
    long m_nIndex, m_nLowerBound, m_nUpperBound;
};


class BadBufferException : public BaseException {
protected:
    void PrintErrMsg() { 
        printf(m_szErrMsg, m_nBufSize, m_nCorrectSize); };
        
public:
    BadBufferException(unsigned long nBufSize, unsigned long nCorrectSize) {
        m_nBufSize = nBufSize;
        m_nCorrectSize = nCorrectSize;
        m_szErrMsg = "Buffer incorrect size of %d, should be %d";
    };
    
    unsigned long m_nBufSize, m_nCorrectSize;
};


#endif
