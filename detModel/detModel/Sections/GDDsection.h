#ifndef GDDSECTION_H
#define GDDSECTION_H
#include <vector>

#include "detModel/Management/GDDsectionsVisitor.h"

class GDDchoice;
class GDDvolume;
class GDDshape;

/**
 * @author D.Favretto & R.Giannitrapani
 */
class GDDsection {
public:
  GDDsection(std::string pname,
	     std::string pversion,
	     std::string pdate,
	     std::string pauthor,
	     std::string ptopVolumeRef):name(pname),version(pversion),date(pdate),
    author(pauthor),topVolumeRef(ptopVolumeRef),
    topVolume(0){};
  
  
  ~GDDsection();

  void addVolume(GDDvolume* nextVolume){volumes.push_back(nextVolume);}
  void addChoice(GDDchoice* nextChoice){choices.push_back(nextChoice);}
  void Accept(GDDsectionsVisitor* v);
  void AcceptNotRec(GDDsectionsVisitor* v);

  std::string getName(){return name;}    
  std::string getVersion()const{return version;};
  std::string getDate()const{return date;};
  std::string getAuthor()const{return author;};
  GDDvolume* getTopVolume()const{return topVolume;};
  std::string getTopVolumeRef()const{return topVolumeRef;};

  void setTopVolume(GDDvolume* pvolume){topVolume = pvolume;};
  void setVersion(std::string pver){version = pver;};
  void setDate(std::string pdate){date = pdate;};
  void setAuthor(std::string pauth){author = pauth;};
  void setTopVolumeRef(std::string pref){topVolumeRef = pref;}

  std::vector< GDDvolume * > getVolumes()const{return volumes;}
  std::vector < GDDchoice * > getChoices()const{return choices;}



 private:
  /** @link aggregation
   *  @supplierCardinality 1..n */
  std::vector < GDDvolume * > volumes;
  /** @link aggregation */
  std::vector < GDDchoice * > choices;
  std::string name;
  std::string version;
  std::string date;
  std::string author;
  std::string DTDVersion;
  std::string topVolumeRef;
  GDDvolume* topVolume;
};
#endif //GDDSECTION_H

