#ifndef GDDSECTION_H
#define GDDSECTION_H
#include <vector>
#include "detModel/Sections/GDDvolume.h"
#include "detModel/Management/GDDsectionsVisitor.h"

class GDDchoice;

enum logicaltpe{unions,intersection,subtraction};

/**
 * @author D.Favretto & R.Giannitrapani
 */
class GDDsection {
public:
  GDDsection(string pname,
	     string pversion,
	     string pdate,
	     string pauthor,
	     string ptopVolumeRef):name(pname),version(pversion),date(pdate),
    author(pauthor),topVolumeRef(ptopVolumeRef),
    topVolume(0){};
  
  virtual ~GDDsection();
  /**
   *
   */ 
  void addVolume(GDDvolume* nextVolume){volumes.push_back(nextVolume);}
  /**
   *
   */ 
  void addChoice(GDDchoice* nextChoice){choices.push_back(nextChoice);}
  /**
   *
   */ 
  void Accept(GDDsectionsVisitor* v);
  /**
   *
   */ 
  void AcceptNotRec(GDDsectionsVisitor* v);
  /**
   *
   */ 
  string getName(){return name;}    
  /**
   *
   */ 
  string getVersion()const{return version;};
  /**
   *
   */ 
  string getDate()const{return date;};
  /**
   *
   */ 
  string getAuthor()const{return author;};
  /**
   *
   */ 
  GDDvolume* getTopVolume()const{return topVolume;};
  /**
   *
   */ 
  string getTopVolumeRef()const{return topVolumeRef;};
  /**
   *
   */ 
  void setTopVolume(GDDvolume* pvolume){topVolume = pvolume;};
  /**
   *
   */ 
  void setVersion(string pver){version = pver;};
  /**
   *
   */ 
  void setDate(string pdate){date = pdate;};
  /**
   *
   */ 
  void setAuthor(string pauth){author = pauth;};
  /**
   *
   */ 
  void setTopVolumeRef(string pref){topVolumeRef = pref;}
  /**
   *
   */
  vector< GDDvolume * > getVolumes()const{return volumes;}
  /**
   *
   */
  vector < GDDchoice * > getChoices()const{return choices;}
 private:
  /** @link aggregation
   *  @supplierCardinality 1..n */
  vector < GDDvolume * > volumes;
  /** @link aggregation */
  vector < GDDchoice * > choices;
  string name;
  string version;
  string date;
  string author;
  string DTDVersion;
  string topVolumeRef;
  GDDvolume* topVolume;
};
#endif //GDDSECTION_H

