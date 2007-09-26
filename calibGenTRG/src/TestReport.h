// class to write test reports similar to python online reports
// M. Kocian,SLAC, 9/20/06
//
#include <fstream>

class TestReport{
 public:
  TestReport(char* filename);
  void newheadline(char*);
  void additem(char*,char*);
  void addlink(char*,char*,char*);
  void addstatus(char*);
  void starttable( char*[],int);
  void addtableline(char*[],int);
  void endtable();
  void writereport();  
  void greentext(char* , const char*);
  void redtext(char* , const char*);
  void linktext(char*, const char*, const char*);
 private:
  std::ofstream of;
};
