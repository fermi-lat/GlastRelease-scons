// $Header$
#ifndef MOOT_DBUTILEXCEPTION_H
#define MOOT_DBUTILEXCEPTION_H
#include <exception>

namespace MOOT {

  class DbUtilException : std::exception {
  public:
    DbUtilException(const std::string& extraInfo = "", int code=0) : 
      std::exception(),
      m_name("DbUtilException"), m_extra(extraInfo), m_code(code) {}
    virtual ~DbUtilException() throw() {}
    virtual std::string getMsg() {
      std::string msg = m_name + ": " + m_extra;
      return msg;}
    virtual int getCode() const { return m_code;}
    virtual const char* what() {
      return m_extra.c_str();
    }
  protected: 
    std::string m_name;
  private:
    std::string m_extra;
    int         m_code;
  };

}
#endif
