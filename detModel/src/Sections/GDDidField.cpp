#include <string>
#include "detModel/Sections/GDDidField.h"


void GDDidField::setName(string pname){
  if (pname=="layer")name=layer;
  else if (pname=="tower")name=tower;
  else if (pname=="tray")name=tray;
  else if (pname=="row")name=row;
  else if (pname=="col")name=col;
  else if (pname=="topside")name=topside;
  else if (pname=="xtal")name=xtal;
  else if (pname=="orient")name=orient;
  else if (pname=="nsew")name=nsew;
  else if (pname=="diode")name=diode;
  else if (pname=="wafer")name=wafer;
  else if (pname=="SiPosition")name=SiPosition;
}

string GDDidField::getName(){
  std::string s;
  if (name==layer)s="layer";
  else if (name==tower)s="tower";
  else if (name==tray)s="tray";
  else if (name==row)s="row";
  else if (name==col)s="col";
  else if (name==topside)s="topside";
  else if (name==xtal)s="xtal";
  else if (name==orient)s="orient";
  else if (name==nsew)s="nsew";
  else if (name==diode)s="diode";
  else if (name==wafer)s="wafer";
  else if (name==SiPosition)s="SiPosition";
  return s;
}
