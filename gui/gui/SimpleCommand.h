// $Id$
// template for simple commands

#ifndef CONTROL_SIMPLECOMMAND_H
#define CONTROL_SIMPLECOMMAND_H

#include "gui/Command.h"

namespace gui {

template <class Receiver>
class SimpleCommand : public Command {
	// A template class for simple commands. The template parameter is the
	// class which is the receiver, i.e the class which actually performs
	// the operations.
public:
   typedef void (Receiver::* Action)();

   SimpleCommand(Receiver* r, Action a):m_receiver(r),m_action(a){};
   // Constructor needs the receiver class and the method to invoke in that class

   virtual void execute();
   // All a command can do - tells the receiver to perform the action.

private:
   Receiver* m_receiver;
   Action m_action;

};

template <class Receiver> inline
void SimpleCommand<Receiver>::execute(){
	(m_receiver->*m_action)();
}

} // namespace gui

#endif // CONTROL_SIMPLECOMMAND_H
