
#include "LogText.h"

LogText::LogText(FXComposite *p,FXObject* tgt,FXSelector sel,FXuint opts,FXint x,FXint y,FXint w,FXint h):
  FXText(p, tgt, sel, opts, x, y, w, h)
{
  // Nothing else to do other than calling the parent constructor
}


// Method for using FXText as a log window
void LogText::logText(FXString text)
{
  appendText(text.text(), text.length());
  makePositionVisible(lineStart(length));
}
