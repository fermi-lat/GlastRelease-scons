
#ifndef FileHeader_H
#define FileHeader_H

#include <TObject.h>
#include <TString.h>
#include <TMap.h>
class ISvcLocator ;

/** @class FileHeader
 * @brief Set of dynamically named attributes, storable inside a root file.
 *
 * $Header$
 */

class FileHeader: public TObject {

public:

    FileHeader() ;
    virtual ~FileHeader() ;

    void DeleteAll();
 
    void setInteger( const char * name, Int_t ) ;
    Int_t getInteger( const char * name ) const ;

    void setDouble( const char * name, Double_t ) ;
    Double_t getDouble( const char * name ) const ;

    void setString( const char * name, const TString & ) ;
    TString getString( const char * name ) const ;

private:

    TMap m_integers ;
    TMap m_doubles ;
    TMap m_strings ;

    ClassDef(FileHeader,1)
}; 
 
#endif





