/*+---------------------- Copyright notice ----------------------------+*/
/*| Copyright (C) 1995, Guy Barrand, LAL Orsay, (barrand@lal.in2p3.fr) |*/
/*|  Permission to use, copy, modify, and distribute this software     |*/
/*| and its documentation for any purpose and without fee is hereby    |*/
/*| granted, provided that the above copyright notice appear in all    |*/
/*| copies and that both that copyright notice and this permission     |*/
/*| notice appear in supporting documentation.  This software is       |*/
/*| provided "as is" without express or implied warranty.              |*/
/*+---------------------- Copyright notice ----------------------------+*/

#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "CalRecon/Midnight.h"

#define changeBlockSize realloc
#define freeBlock free
#define allocateBlock malloc

static void DeleteString(char*);
static char* DuplicateString(const char*);
static char* ConcatenateString(char*,const char*);
static char* Create(int);

#define MINIMUM(a,b) ((a)<(b)?a:b)

//////////////////////////////////////////////////////////////////////////////
MString::MString(
)
:fString(NULL)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  fString = DuplicateString("");
}
//////////////////////////////////////////////////////////////////////////////
MString::MString(
 const char* aString
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  fString = DuplicateString(aString==NULL?"":aString);
}
//////////////////////////////////////////////////////////////////////////////
MString::MString(
 const MString& aFrom
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  fString = DuplicateString(aFrom.fString);
}
//////////////////////////////////////////////////////////////////////////////
MString::~MString(
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  DeleteString(fString);
}
//////////////////////////////////////////////////////////////////////////////
const char* MString::data(
) const
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  return fString;
}
//////////////////////////////////////////////////////////////////////////////
int MString::length(
) const
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  if(fString==NULL) return 0;
  return strlen(fString);
}
//////////////////////////////////////////////////////////////////////////////
void MString::resize(
 int aLength
)
//////////////////////////////////////////////////////////////////////////////
// Truncate or add blanks as necessary.
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  if(fString==NULL) return;
  if(aLength<0) return;
  int l = strlen(fString);
  if(aLength<l) {
    fString[aLength] = '\0';
  } else {
    char* s = Create(aLength-l);
    if(s!=NULL) {
      fString = ConcatenateString(fString,s);
      freeBlock(s);
    }
  }
}
//////////////////////////////////////////////////////////////////////////////
MString& MString::replace(
 int aStart
,int aLength
,const MString& aString
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  if(fString==NULL) return *this;
  if(aString.fString==NULL) return *this;
  if( (aStart<0) || (aStart>=length()) ) return *this;
  int begin = aStart;
  int end = MINIMUM(aStart + aLength - 1,length()-1);
  int l = aString.length();
  int pos = 0;
  for(int count=begin;count<=end;count++,pos++) {
    if(pos<l) 
      fString[count] = aString.fString[pos];
    else
      fString[count] = ' ';
  }
  return *this;
}
//////////////////////////////////////////////////////////////////////////////
MString& MString::operator +=(
 const char* aString
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  fString = ConcatenateString(fString,aString);
  return *this;  
}
//////////////////////////////////////////////////////////////////////////////
MString& MString::operator +=(
 const MString& aString
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  *this += aString.fString;
  return *this;
}
//////////////////////////////////////////////////////////////////////////////
MString& MString::operator =(
 const char* aString
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  if((fString!=NULL) && (aString==fString) ) return *this; 
  DeleteString(fString);
  fString = DuplicateString(aString==NULL?"":aString);
  return *this;
}
//////////////////////////////////////////////////////////////////////////////
MString& MString::operator =(
 const MString& aFrom
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  return (*this = aFrom.fString);
}
//////////////////////////////////////////////////////////////////////////////
MString& MString::operator =(
 char aChar
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  DeleteString(fString);
  fString = Create(1);
  fString[0] = aChar;
  return *this;
}
//////////////////////////////////////////////////////////////////////////////
MString::operator const char*(
) const 
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{ 
  return fString;
}
//////////////////////////////////////////////////////////////////////////////
char& MString::operator[](
 int aIndex
)
//////////////////////////////////////////////////////////////////////////////
// To do x[i] = c;
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{ 
  if((aIndex<0)||(aIndex>=length())) {
    printf("MString::operator[] : bad index %d %d\n",aIndex,length());
    //exit(EXIT_FAILURE);
    aIndex = 0;
  }
  return fString[aIndex]; 
}
//////////////////////////////////////////////////////////////////////////////
char MString::operator[](
 int aIndex
) const 
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{ 
  if((aIndex<0)||(aIndex>=length())) {
    printf("MString::operator[] : bad index %d %d\n",aIndex,length());
    //exit(EXIT_FAILURE);
    aIndex = 0;
  }
  return fString[aIndex]; 
}
//////////////////////////////////////////////////////////////////////////////
char& MString::operator()(
 int aIndex
)
//////////////////////////////////////////////////////////////////////////////
// To do x(i) = c;
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{ 
  return (*this)[aIndex]; 
}
//////////////////////////////////////////////////////////////////////////////
char MString::operator()(
 int aIndex
) const 
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{ 
  return (*this)[aIndex]; 
}
//////////////////////////////////////////////////////////////////////////////
MString MString::operator()(
 int aStart
,int aLength
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  MString s;
  if(fString==NULL) return s;
  if( (aStart<0) || (aStart>=length()) ) return s;
  int begin = aStart;
  int end = MINIMUM(aStart + aLength - 1,length()-1);
    DeleteString(s.fString);
  s.fString = Create(end - begin + 1);
  int pos = 0;
  for(int count=begin;count<=end;count++,pos++) {
    s.fString[pos] = fString[count];
  }
  return s;
}
//////////////////////////////////////////////////////////////////////////////
int operator ==(
 const MString& a1
,const MString& a2
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  if( a1.fString ==  a2.fString) return 1; 
  if( (a1.fString==NULL) || (a2.fString==NULL) ) return 0; 
  return (strcmp(a1.fString,a2.fString)==0 ? 1: 0); 
}
//////////////////////////////////////////////////////////////////////////////
int operator !=(
 const MString& a1
,const MString& a2
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  return (a1==a2 ? 0 : 1);
}
//////////////////////////////////////////////////////////////////////////////
int operator ==(
 const MString& a1
,const char* a2
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  if( a1.fString ==  a2) return 1; 
  if( (a1.fString==NULL) || (a2==NULL) ) return 0; 
  return (strcmp(a1.fString,a2)==0 ? 1: 0); 
}
//////////////////////////////////////////////////////////////////////////////
int operator !=(
 const MString& a1
,const char* a2
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  return (a1==a2 ? 0 : 1);
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
char* ConcatenateString (
 char* aString
,const char* aFrom 
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  if(aFrom==NULL) return aString;
  if(aString==NULL) {
    aString = DuplicateString(aFrom);
  } else if(*aString=='\0') {
    DeleteString(aString);
    aString = DuplicateString(aFrom);
  } else { 
    int lto = strlen(aString);
    int lfrom = strlen(aFrom);
    int length = lto+lfrom;
    char* str = 
      (char*)changeBlockSize(aString,(size_t)((length+1)*sizeof(char)));
    if(str==NULL) return aString;
    str[length] = '\0';
    aString = str;
    strncpy(aString+lto,aFrom,lfrom);
    aString[lto+lfrom] = '\0';
  }
  return aString;
}
//////////////////////////////////////////////////////////////////////////////
static void DeleteString(
 char* This
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  if(This==NULL) return;
  freeBlock(This);
}
//////////////////////////////////////////////////////////////////////////////
static char* DuplicateString(
 const char* This
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  if(This==NULL) return NULL;
  int length = strlen(This);
  char* string = (char*)allocateBlock((length+1)*sizeof(char));
  if(string==NULL) return NULL;
  string[length] = '\0';
  return strcpy(string,This);
}
//////////////////////////////////////////////////////////////////////////////
static char* Create(
 int a_length
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  if(a_length<=0) a_length = 0;
  char* string = (char*)allocateBlock((a_length+1)*sizeof(char));
  if(string==NULL) return NULL;
  string[a_length] = '\0';
  for(int count=0;count<a_length;count++) string[count] = ' ';
  return string;
}

