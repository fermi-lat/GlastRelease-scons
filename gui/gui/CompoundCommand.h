//     $Id$
// Project: Atlas top-level prototype
//  Author: Toby Burnett
//

#ifndef COMPOUNDCOMMAND_H
#define COMPOUNDCOMMAND_H

#include "gui/Command.h"
#include <vector>

namespace gui {

typedef  std::vector<Command*  > CommandList;
    
class CompoundCommand : public Command , private CommandList {
    // A Command which is a list of Commands
    
public:
    CompoundCommand(Command* cmd = 0);
    // constructor can insert the first command	
    
    ~CompoundCommand();
    
    void append(Command* cmd);
    // add a command to the list
    
    virtual void execute();
    // Like all commands, it needs an execute method.
};

} //namespace gui

#endif

