//   $Header$
// Project:
//  author: Toby Burnett
//
//  Declare text output window for Motif environment, with
//  possibility of clear and label of window

#ifndef XOSTREAM_H
#define XOSTREAM_H
 

#include "gui/GUI.h"


class Xostream : public std::ostream
{
public:
    Xostream(const char* name);
    void clear();
    void flush();
private:
};


#endif

