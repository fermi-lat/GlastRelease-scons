
#ifndef FileHeader_H
#define FileHeader_H

#include <TObject.h>
#include <TString.h>
#include <TMap.h>
class ISvcLocator ;

/*! 

 @class FileHeader
 @brief Set of dynamically named attributes, storable inside a root file.

 Each ROOT file should contains one instance of this class, with
 data concerning all the file events. The interface is very open: anyone
 can store there any string, integer or double associated to
 a name.
 
 @author David Chamont - CNRS IN2P3 LLR Ecole Polytechnique

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

    TMapIter& getIntMapIter(Bool_t dir=kIterForward) const; 
 
    TMapIter& getDoubleMapIter(Bool_t dir=kIterForward) const ; 

    TMapIter& getStringMapIter(Bool_t dir=kIterForward) const ; 

private:

    TMap m_integers ;
    TMap m_doubles ;
    TMap m_strings ;

    ClassDef(FileHeader,1)
}; 
 
#endif





