#ifndef SEG_H
#define SEG_H

namespace detModel{

  /**
   * @author D.Favretto & R.Giannitrapani
   */
  class Seg {
  public:
    enum axisDir{xDir,yDir,zDir};

  public:
  
    Seg():nSeg(1){;};
    /**
     *This method sets the value of the attibute reason
     */ 
    void setReason(std::string preason){reason=preason;}
    /**
     *This method returns the value of the attribute reason
     */ 
    std::string getReason()const{return reason;}
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
    void setAxis(std::string paxis){if (paxis=="X")axis=xDir;
    else if(paxis=="Y")axis=yDir;
    else axis=zDir;
    }
    /**
     *This method returns the value of the attribute axis
     */ 
    axisDir getAxis()const{return axis;}
  
  private:    
    axisDir axis;
    std::string reason;
    int nSeg; ///default value 1
  };
}
#endif //SEG_H
