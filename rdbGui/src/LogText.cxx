
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

std::ostringstream* LogText::getOutStream()
{
  return &m_outStream;
}

std::ostringstream* LogText::getErrStream()
{
  return &m_errStream;
}

void LogText::update()
{
  logText(m_errStream.str().c_str());
  logText(m_outStream.str().c_str());
  m_errStream.flush();
  m_outStream.flush();
}
