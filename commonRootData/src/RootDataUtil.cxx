#include "commonRootData/RootDataUtil.h"

void rootdatautil::CheckBounds(long nIndex, long nLower, long nUpper) {
	if ((nIndex < nLower) || (nIndex > nUpper))
		throw(IndexOutOfBoundsException(nIndex, nLower, nUpper));
}

void rootdatautil::CheckBufSize(unsigned long nBufSize, unsigned long nCorrectSize) {
	if (nBufSize != nCorrectSize)
		throw(BadBufferException(nBufSize, nCorrectSize));
}
