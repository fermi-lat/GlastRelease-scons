#ifndef GDDSEG_H
#define GDDSEG_H
enum axisType{X,Y,Z};

/**
 * @author D.Favretto & R.Giannitrapani
 */
class GDDseg {
 public:
  
  GDDseg():nSeg(1){;};
  /**
   *This method sets the value of the attibute reason
   */ 
  void setReason(string preason){reason=preason;}
  /**
   *This method returns the value of the attribute reason
   */ 
  string getReason()const{return reason;}
  /**
   *This method sets the value of the attibute nSeg
   */ 
  void setnSeg(int pnSeg){nSeg=pnSeg;}
  /**
   *This method returns the the value of the attribute nSeg
   */ 
  int getnSeg()const{return nSeg;}
  /**
   *This method sets the value of the attribute axis 
   */ 
  void setAxis(string paxis){if (paxis=="X")axis=X;
  else if(paxis=="Y")axis=Y;
  else axis=Z;
  }
  /**
   *This method returns the value of the attribute axis
   */ 
  axisType getAxis()const{return axis;}
  
  private:    
  axisType axis;
  string reason;
  int nSeg; ///default value 1
  };
#endif //GDDSEG_H
  
  
  
