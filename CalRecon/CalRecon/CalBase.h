
#ifndef CALBASE_H
#define CALBASE_H

//############################
class CalBase
//############################
{
public:
	enum SIDE  {NEG, POS};
	enum RANGE {LEX,LE,HEX,HE};

	static SIDE  side(int i);
	static RANGE range(int i);
};


const int CALNRANGES = 4;
const int CALNSIDES  = 2;

#endif
