//     $Id$
//  Author: Toby Burnett


#ifndef LOOPCOMMAND_H
#define LOOPCOMMAND_H

#include "gui/Command.h"
namespace gui {

class LoopCommand : public Command {
    // Loops over a contained command.
    // It provides functions to control the loop, and will return 
    // <Command>> instances for each.
public:
    LoopCommand(Command* command, long count);
    // Constructor has a <<Command>> to loop over for "count" times

    ~LoopCommand();
    // dtor - delete the command

    void setCount(long count){m_count = count;}
    // Set the loop count. It will have been set initially by the 
    // constructor.

    int getCount()const{return m_count;}
    // Return the loop count (includes those executed and those yet
    // to be done)

    int getCurrent()const{return m_current;}
    // Return the number of loops still to be executed

    void setCommand(Command* cmd){m_command =cmd;}
    // Set a new command

    void execute();
    // Start looping for the specified number of times.

    void pause();
    // Pause the loop
    
    void abort();
    // Abort the loop

    void resume();
    // Resume looping

    void stop();
    // Stop the loop  
 
    enum Status {STOPPED, RUNNING, PAUSED, ABORTED};

    Status status()const;
    // Access to the status. It may be STOPPED, RUNNING, PAUSED or ABORTED

    Command* abortCommand();
    // Return the command that will abort the loop
    
    Command* pauseCommand();
    // Return the command that will pause the loop

    Command* resumeCommand();
    // Return the command that will resume the loop

private:
    class LCommand : public Command {
	// inner private class that defines simple commands affecting the loop itself  
	friend class LoopCommand;
	typedef void (LoopCommand::* Action)();

	LCommand(LoopCommand* r, Action a):
	    _receiver(r), _action(a) {}

	virtual void execute(){ (_receiver->*_action)();}

	LoopCommand* _receiver;
	Action _action;
    };

    Command*  m_command; // the command to execute in the loop
    unsigned long      m_count;   // max loop count, set by constructor
    unsigned long      m_current; // loop counter: count up to the max
    Status    m_status;  // control start/stop
};
} //namespace dove
#endif
