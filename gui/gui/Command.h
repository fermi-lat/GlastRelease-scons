// Command.h,v 1.1 1995/05/15 17:39:41 burnett Exp
// Abstract base class for commands

#ifndef CONTROL_COMMAND_H
#define CONTROL_COMMAND_H

namespace gui {

class Command {
	// This is based on the Command pattern from "Design Patterns".
	// Command has just a pure virtual function execute() for which
	// concrete classes must provide an implementation.

public:
    virtual ~Command();
	// virtual destructor for a base class

	virtual void execute() = 0;
    // subclasses must implement execute

protected:
    Command();
};

inline 
Command::~Command() {}

inline 
Command::Command() {}

} //namespace gui
#endif // CONTROL_COMMAND_H

