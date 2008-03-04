// class to write test reports similar to python online reports
// M. Kocian,SLAC, 9/20/06
//

#include <iostream>
#include "TestReport.h"
#include <string>

TestReport::TestReport(char *filename){
  of.open(filename);  
  of<<"<HTML>"<<std::endl;
of<<"<HEAD>"<<std::endl;
of<<"<META http-equiv=\"Content-Type\" content=\"text/html; charset=UTF-8\">"<<std::endl;
of<<"<TITLE>TRG Test Report</TITLE>"<<std::endl;
of<<"<STYLE media=\"print\" type=\"text/css\">#toc {display:none}p {font: 12pt serif; page-break-inside:avoid}h1,h2 {page-break-before:always}</STYLE>"<<std::endl;
of<<"<STYLE media=\"screen\" type=\"text/css\">h1,h2 {border-top-style:solid}</STYLE>"<<std::endl;
of<<"</HEAD>"<<std::endl;
of<<"<BODY>"<<std::endl;
of<<"<Title>TRG Test Report</Title>"<<std::endl;
of<<"<cdatanode>"<<std::endl;
of<<"    <STYLE TYPE=\"text/css\">"<<std::endl;
of<<"    H1 { font-size: x-large; color: blue }"<<std::endl;
of<<"    H2 { font-size: large; color: blue }"<<std::endl;
of<<"    H3 { font-size: medium; color: blue }"<<std::endl;
of<<"    </STYLE>"<<std::endl;
of<<""<<std::endl;
of<<"    </cdatanode>"<<std::endl;
of<<"<p>&nbsp;</p>"<<std::endl;
of<<"<p>&nbsp;</p>"<<std::endl;
of<<"<p>&nbsp;</p>"<<std::endl;
of<<"<p>&nbsp;</p>"<<std::endl;
of<<"<p>&nbsp;</p>"<<std::endl;
of<<"<p>&nbsp;</p>"<<std::endl;
of<<"<p>&nbsp;</p>"<<std::endl;
of<<"<center>"<<std::endl;
of<<"<b><font size=\"6\">"<<std::endl;
of<<"<Line no=\"1\">"<<std::endl;
of<<"      <cdatanode>GLAST LAT Trigger</cdatanode>"<<std::endl;
of<<"    </Line>"<<std::endl;
of<<""<<std::endl;
of<<"</font></b>"<<std::endl;
of<<"</center>"<<std::endl;
of<<"<center>"<<std::endl;
of<<"<b><font size=\"6\">"<<std::endl;
of<<"<Line no=\"2\">"<<std::endl;
of<<"      <cdatanode>Test Report</cdatanode>"<<std::endl;
of<<"    </Line>"<<std::endl;
of<<"</font></b>"<<std::endl;
of<<"</center>"<<std::endl;
of<<"<p>&nbsp;</p>"<<std::endl;
of<<"<p>&nbsp;</p>"<<std::endl;
of<<"<p>&nbsp;</p>"<<std::endl;
of<<"<p>&nbsp;</p>"<<std::endl;
of<<"<p>&nbsp;</p>";

}

void TestReport::newheadline(char* text){
  of<<"<h1> <Caption>"<<text<<"</Caption> </h1>"<<std::endl;
}
void TestReport::additem(char* text1, char* text2){
  of<<"<p> <b><Label>"<<text1<<"</Label>:&nbsp;</b> <Text>"<<text2<<"</Text> </p>"<<std::endl;
}
void TestReport::addlink(char* text1,char* link, char* text2){
  of<<"<p> <b><Label>"<<text1<<"</Label>:&nbsp;</b><a href=\""<<link<<"\"> <Text>"<<text2<<"</Text> </a> </p>"<<std::endl;
}

void TestReport::addstatus(char* text){
  of<<"<p> <b><Label>Status</Label>:&nbsp;</b> <Text> <b>"<<std::endl;
  std::string s(text);
  if (s=="PASSED" || s=="passed" || s=="Passed")of<<"<font color=\"#33CC33\">"<<std::endl;
  else if (s=="FAILED" || s=="failed" || s=="Failed")of<<"<font color=\"#FF0000\">"<<std::endl;
  else if (s=="ABORTED" || s=="aborted" || s=="Aborted")of<<"<font color=\"#FFAA33\">"<<std::endl;
  else of<<"<font color=\"#000000\">"<<std::endl;
  of<<text;
  of<<"</font> </b> </Text> </p>"<<std::endl;
}

void TestReport::starttable(char* headers[],int columns){
  of<<"<Table width=\"100%\" border=\"1\"> <TR>"<<std::endl;
  for (int i=0;i<columns;i++){
    of<<"<TH align=\"center\"> <cdatanode>"<<headers[i]<<"</cdatanode> </TH>"<<std::endl;
  }
  of<<"<TR>"<<std::endl;
}
void TestReport::addtableline(char* line[],int columns){
  of << "<TR>";
  for (int i=0;i<columns;i++){
    of<<"<TD width=\"center\" align=\"center\"> <cdatanode>"<<line[i]<<"</cdatanode> </TD>"<<std::endl;
  }
}

void TestReport::endtable(){
  of<<"</Table>";
}
void TestReport::addimage(char* location, char* description){
  of<<"<img src="<<location<<" alt="<<description<<">"<<std::endl;
}
    
void TestReport::writereport(){
  of<<"</BODY> </HTML>"<<std::endl;
  of.close();  
}
void TestReport::greentext(char* textout,const char* textin){
  sprintf(textout,"<font color=\"#33CC33\">%s</font>",textin);
}
void TestReport::redtext(char* textout,const char* textin){
  sprintf(textout,"<font color=\"#FF0000\">%s</font>",textin);
}
void TestReport::linktext(char * textout,const char* textin,const char* linkin){
  sprintf(textout,"<a href=%s> <Text>%s</Text> </a>",linkin,textin);
}
