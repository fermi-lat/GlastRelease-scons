/*+---------------------- Copyright notice ----------------------------+*/
/*| Copyright (C) 1995, Guy Barrand, LAL Orsay, (barrand@lal.in2p3.fr) |*/
/*|  Permission to use, copy, modify, and distribute this software     |*/
/*| and its documentation for any purpose and without fee is hereby    |*/
/*| granted, provided that the above copyright notice appear in all    |*/
/*| copies and that both that copyright notice and this permission     |*/
/*| notice appear in supporting documentation.  This software is       |*/
/*| provided "as is" without express or implied warranty.              |*/
/*+---------------------- Copyright notice ----------------------------+*/

// Implemented API is STL compatible.

#ifndef MString_h
#define MString_h

class MString {
public:
  MString();
  MString(const char*);
  MString(const MString&);
  ~MString();
  const char* data() const;
  int length() const;
  void resize(int);
  MString& replace(int,int,const MString&);
  MString& operator +=(const char*);
  MString& operator +=(const MString&);
  MString& operator =(const char*);
  MString& operator =(const MString&);
  MString& operator =(char);
  operator const char*() const;
  char& operator[](int);
  char  operator[](int) const;
  char& operator()(int);
  char  operator()(int) const;
  MString operator()(int,int);
  friend int operator ==(const MString&,const MString&);
  friend int operator !=(const MString&,const MString&);
  friend int operator ==(const MString&,const char*);
  friend int operator !=(const MString&,const char*);
private:
  char* fString;
};

#endif




