//     $Id$
//  Author: Toby Burnett
//
// Implementation of CompoundCommand

#include "gui/CompoundCommand.h"
#include <algorithm>

namespace gui {

CompoundCommand::CompoundCommand(Command* cmd)
{
    if(cmd ) append(cmd);
}

void CompoundCommand::append(Command * cmd)
{
    push_back(cmd);
}

CompoundCommand::~CompoundCommand()
{
    iterator ic = begin();
    while (ic != end() )
        delete *ic++;
}

//static void Executor (Command *p) { p->execute (); }
class Executor { 
public: 
    void operator()(Command* c){c->execute();} 
};

void CompoundCommand::execute()
{
    std::for_each (begin (), end (), Executor());
}

} // namespace gui
