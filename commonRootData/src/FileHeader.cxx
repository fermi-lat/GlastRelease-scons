
#include "commonRootData/FileHeader.h"
#include "commonRootData/RootObj.h"
#include <TObjString.h>
#include <TNamed.h>
#include <iostream>

ClassImp(FileHeader)

FileHeader::FileHeader() {
}

FileHeader::~FileHeader() {
	DeleteAll() ;
}

void FileHeader::DeleteAll() {
	m_integers.DeleteAll() ;
	m_doubles.DeleteAll() ;
	m_strings.DeleteAll() ;
}

void FileHeader::setInteger( const char * name, Int_t value ) {
    m_integers.Add(new TNamed(name,""),new RootObj<Int_t>(value)) ;
}

Int_t FileHeader::getInteger( const char * name ) const {
	TObject * object = m_integers.GetValue(name) ;
	if (!object)
	  return -1 ;
	else
	  return ((RootObj<Int_t> *)object)->value() ;
}
 
void FileHeader::setDouble( const char * name, Double_t value ) {
    m_doubles.Add(new TNamed(name,""),new RootObj<Double_t>(value)) ;
}

Double_t FileHeader::getDouble( const char * name ) const {
	TObject * object = m_doubles.GetValue(name) ;
	if (!object)
	  return 0 ;
	else
	  return ((RootObj<Double_t> *)object)->value() ;
}
 
void FileHeader::setString( const char * name, const TString & value ) {
    m_strings.Add(new TNamed(name,""),new TObjString(value)) ;
}

TString FileHeader::getString( const char * name ) const {
	TObject * value = m_strings.GetValue(name) ;
	if (!value)
	  return "Empty String" ;
	else
	  return ((TObjString *)value)->GetString() ;
}

TMapIter& FileHeader::getIntMapIter(Bool_t dir) const {
    static TMapIter m(&m_integers, dir);
    return m;
}

TMapIter& FileHeader::getDoubleMapIter(Bool_t dir) const {
    static TMapIter m(&m_doubles, dir);
    return m;
}

TMapIter& FileHeader::getStringMapIter(Bool_t dir) const {
    static TMapIter m(&m_strings, dir);
    return m;
}
 
