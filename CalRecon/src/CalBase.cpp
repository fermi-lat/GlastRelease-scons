
#include "CalRecon/CalBase.h"

//#####################################
enum CalBase::RANGE CalBase::range(int i)
//#####################################
{
	RANGE r = LEX;
	r = (i == 0? LEX : LE);
	if (i >=2) r = ( i == 2 ? HEX : HE);
	return r;
}
//#####################################
enum CalBase::SIDE CalBase::side(int i)
//#####################################
{
	return ( i == 0? NEG : POS);
}
