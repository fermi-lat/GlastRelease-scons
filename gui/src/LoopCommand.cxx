//     $Id$
//  Author: Toby Burnett
//
// Implementation of the LoopCommand class


#include "gui/LoopCommand.h"
#include <limits.h>

namespace gui {

LoopCommand::LoopCommand(Command* command, long count)
: m_command(command)
, m_count(count)
{
    stop();
}

LoopCommand::~LoopCommand()
{
    delete m_command;
}

 
void LoopCommand::execute()
{
    m_current = 0;
    m_status = PAUSED;
    resume();
}


void LoopCommand::resume()
{
    if( m_status != PAUSED)
	return;
    m_status = RUNNING;
    for(; m_status==RUNNING && m_current< m_count; ++m_current)  {
	  m_command->execute();
    }
    if( m_status==RUNNING)
	  m_status= STOPPED;
}

void LoopCommand::pause()
{
    m_status = PAUSED;
}
void LoopCommand::stop() {   m_status= STOPPED;  m_current = ULONG_MAX;}
void LoopCommand::abort(){   m_status = ABORTED; m_current = ULONG_MAX;}


LoopCommand::Status LoopCommand::status() const
{
	return m_status;
}
Command* LoopCommand::pauseCommand() { return new LCommand(this, &LoopCommand::pause);}
Command* LoopCommand::abortCommand() { return new LCommand(this, &LoopCommand::abort);}
Command* LoopCommand::resumeCommand(){ return new LCommand(this, &LoopCommand::resume);}
	
} // namespace gui
