// SpectrumFactoryTable.h: interface for the SpectrumFactoryTable class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_SPECTRUMFACTORYTABLE_H__E3EDA893_F360_4B4D_9537_0831B21DD02F__INCLUDED_)
#define AFX_SPECTRUMFACTORYTABLE_H__E3EDA893_F360_4B4D_9537_0831B21DD02F__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <map>
#include <string>
#include <map>
#include <list>
#include "FluxSvc/ISpectrumFactory.h"

//! singleton table holding references to all the spectrum classes.  
//! serves to instantiate new classes based on a name.
class SpectrumFactoryTable : public std::map<std::string, const ISpectrumFactory* > 
{
public:
    void addFactory(std::string name, const ISpectrumFactory* factory ) {
        insert(std::make_pair<std::string, const ISpectrumFactory*>(name, factory));
    }
    
    //! get a new Spectrum object by name with optional parameter list
    Spectrum* instantiate(const std::string& name, const std::string& params) const ;
    Spectrum* instantiate(const std::string& name) const ;
    
    static SpectrumFactoryTable* instance(){
        return (s_instance==0)? new SpectrumFactoryTable : s_instance;
    }
    
    
    std::list<std::string> spectrumList()const;
    
private:
    static SpectrumFactoryTable* s_instance;
    SpectrumFactoryTable(){s_instance=this;}
    
        
};


#endif // !defined(AFX_SPECTRUMFACTORYTABLE_H__E3EDA893_F360_4B4D_9537_0831B21DD02F__INCLUDED_)
