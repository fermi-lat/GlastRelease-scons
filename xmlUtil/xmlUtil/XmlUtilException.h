// $Header$
#ifndef XMLUTIL_XMLUTILEXCEPTION_H 
#define XMLUTIL_XMLUTILEXCEPTION_H
#include <string>
#include <exception>

namespace xmlUtil {
  //! Base exception class for Arith exceptions
  class XmlUtilException : public std::exception {
  public:
    XmlUtilException() : std::exception() {
      m_msg = std::string("XmlUtilException");
      m_name=std::string("XmlUtilException");
    }
    
    virtual ~XmlUtilException() throw() {}
    virtual std::string getMsg() {return m_msg;}
    virtual const char * what() {return m_msg.c_str();}
  private:
    std::string m_msg;
  protected:
    std::string m_name;
    void makeMsg(std::string reason) {m_msg = m_name + ":  " + reason;}
  };
}


#endif
