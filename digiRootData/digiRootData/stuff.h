#ifndef stuff_h
#define stuff_h

////////////////////////////////
//
// stuff.h
//
// Author: Daniel Flath 
// Date:   2000/03/01
//
////////////////////////////////

namespace stuff {

////////////////////////////////////////////////////////////////
//
// Function: Round
//
// Input:
//   float num:  input number to be rounded
//
// Return value:
//   The rounded integer equivalent of num.
//
////////////////////////////////////////////////////////////////
int Round(float num) {
  int retVal;
  int bNegative;

  if (num < 0) {
    num *= -1;
    bNegative = 1;
  } else {
    bNegative = 0;
  }

  int numFloor = (int)num;
  int numCeil = (int)(num + 1);
  if ((num - (float)numFloor) < ((float)numCeil - num))
    retVal = numFloor;
  else
    retVal =  numCeil;

  return bNegative ? -retVal : retVal;
}


////////////////////////////////////////////////////////////////
//
// Function: IntToStr
//
// Input:
//   int in:     number to convert
//   char *buf:  string to store result in
//   int length: size of buf
//
// Return value:
//   Number of characters written to buf excluding the null
//   terminator.
//
////////////////////////////////////////////////////////////////
int IntToStr(int in, char *buf, int length) {
    if (length > 1) {
        // simple case:
        if (in == 0) {
            buf[0] = '0';
            buf[1] = 0;
            return 1;
        }

        int len = 0;
        int temp = in;
        if (temp == 0)
            len = 1;
        while (temp != 0) {
            temp /= 10;
            len++;
        }
        if (in < 0) {
            len++;
        }

        len = (len < length - 1) ? len : length - 1;

        char *ch = buf + len - 1;
        
        temp = in;
        if (temp < 0)
            temp *= -1;

        int chint = 0;
        while (temp != 0) {
            chint = temp % 10;
            temp /= 10;
            ch[0] = chint + 48;
            ch--;
        } 
        
        buf[len] = 0;
        if (in < 0) {
            buf[0] = '-';            
        }

        return len;
    }
    else
        return 0;
}

} // end namespace stuff
#endif
