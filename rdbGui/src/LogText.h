
#ifndef LOGTEXT_H
#define LOGTEXT_H

#include "fx.h"


class LogText : public FXText 
{

 public:

  LogText(FXComposite *p,FXObject* tgt=NULL,FXSelector sel=0,FXuint opts=0,FXint x=0,FXint y=0,FXint w=0,FXint h=0);

  /// Method for using FXText as a log window
  void  logText(FXString text);

};


#endif // end of LOGTEXT_H
