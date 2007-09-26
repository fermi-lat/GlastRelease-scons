#include "eval.h"

using std::map;
using std::vector;
using std::string;


const vector<string> eval::readlistortuple(const string l){
  vector<string> m;
  unsigned int p=0;
  string entry;
  int nquotes, nparen, nbracket, nbrace;
  nquotes=nparen=nbracket=nbrace=0;
  int beg,end;
  char endchar;
  beg=end=endchar=0;
  for (unsigned int i=0;i<l.length();i++){
    if (l[i]=='('){ 
      p=i+1;
      endchar=')'; // tuple
      break;
    }
    if ( l[i]=='['){
      p=i+1;
      endchar=']'; //list
      break;
    }
    if (l[i]!=' ')return m;
  }
  while(1){
    beg=p;
    while (l[beg]==' ')beg++;
    while ((l[p]!=',' && l[p]!=endchar) || nquotes!=0 || nparen!=0 || nbracket!=0 || nbrace!=0){
      if(l[p]=='(')nparen++;
      if(l[p]==')')nparen--;
      if(l[p]=='[')nbracket++;
      if(l[p]==']')nbracket--;
      if(l[p]=='{')nbrace++;
      if(l[p]=='}')nbrace--;
      if(l[p]=='\''&& nquotes==0)nquotes=1;
      if(l[p]=='\''&& nquotes==1)nquotes=0;
      p++;
      if (p>=l.length())return m; 
    }
    end=p;
    while(l[end]==' ')end--;
    entry=l.substr(beg,end-beg); 
    m.push_back(entry);
    if (l[p]==endchar) return m;
    p++;
  }
}

map<string,string> eval::readdict(const string l){
  map<string, string> m;
  unsigned int p=0;
  //const char* l=sl.c_str();
  p=find(l,'{',p);
  if (p==(unsigned)-1)return m;
  int nquotes, nparen, nbracket, nbrace;
  nquotes=nparen=nbracket=nbrace=0;
  int beg, end;
  string key;
  string val;
  while(1){
  beg=p;
  while (l[beg]==' ')beg++;
  while (l[p]!=':' || nquotes!=0 || nparen!=0 || nbracket!=0 || nbrace!=0){
    if(l[p]=='(')nparen++;
    if(l[p]==')')nparen--;
    if(l[p]=='[')nbracket++;
    if(l[p]==']')nbracket--;
    if(l[p]=='{')nbrace++;
    if(l[p]=='}')nbrace--;
    if(l[p]=='\''&& nquotes==0)nquotes=1;
    if(l[p]=='\''&& nquotes==1)nquotes=0;
    p++;
    if (p>=l.length())return m; 
  }
  end=p;
  while(l[end]==' ')end--;
  key=l.substr(beg,end-beg);
  //cout<<"end-beg"<<end-beg<<endl;
  //  cout<<key<<endl;
  p++;
  beg=p;
  while (l[beg]==' ')beg++;
  //cout<<nquotes<<" "<<nparen<<" "<<nbracket<<" "<<nbrace<<endl;
  //cout<<p<<endl;
  while ((l[p]!=',' && l[p]!='}')|| nquotes!=0 || nparen!=0 || nbracket!=0 || nbrace!=0){
    if(l[p]=='(')nparen++;
    if(l[p]==')')nparen--;
    if(l[p]=='[')nbracket++;
    if(l[p]==']')nbracket--;
    if(l[p]=='{')nbrace++;
    if(l[p]=='}')nbrace--;
    if(l[p]=='\''&& nquotes==0)nquotes=1;
    if(l[p]=='\''&& nquotes==1)nquotes=0;
    //      cout<<l[p]<<" "<<nbrace<<" "<<nquotes<<" "<<nparen<<" "<<nbracket<<" "<<nbrace<<endl;
    p++;
    if (p>=l.length()){
      //cout<<nquotes<<" "<<nparen<<" "<<nbracket<<" "<<nbrace<<endl;
      return m; 
    }
  }
  end=p;
  while(l[end]==' ')end--;
  val=l.substr(beg,end-beg);
  //cout<<end-beg<<endl;
  //cout<<val<<endl;
  m[key]=val;
  p++;
}
}
const string eval::unquote(const string l){
  //const char* l=sl.c_str();
  int beg,end;
  beg=end=0;
  for (unsigned int i=0;i<l.length();i++){
    if (l[i]=='\''){
      beg=i+1;
      break;
    }
    if (l[i]!=' '){
      beg=i;
      break;
    }
  }
  for (int i=l.length()-1;i>=0;i--){
    if (l[i]=='\''){
      end=i-1;
      break;
    }
    if (l[i]!=' '){
      end=i;
      break;
    }
  }
  return l.substr(beg,end-beg+1);
  //  nl[end-beg+1]='\0';
}

int eval::find(const string s,char c,int i){
 for (unsigned int j=i;j<s.length();j++){
   if (s[j]==c)return j+1;
 }
 return -1;
}
   
