
#ifndef LOGTEXT_H
#define LOGTEXT_H

#include "fx.h"

#include <sstream>


class LogText : public FXText 
{
 private:
  std::ostringstream m_outStream,m_errStream; 


 public:

  LogText(FXComposite *p,FXObject* tgt=NULL,FXSelector sel=0,FXuint opts=0,FXint x=0,FXint y=0,FXint w=0,FXint h=0);

  /// Method for using FXText as a log window
  void  logText(FXString text);
  
  std::ostringstream* getOutStream();
  std::ostringstream* getErrStream();
  
  void update();

};


#endif // end of LOGTEXT_H
