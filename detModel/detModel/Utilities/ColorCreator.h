#ifndef COLORCREATOR_H
#define COLORCREATOR_H
#include <vector>
#include "detModel/Utilities/Purge.h"

namespace detModel{

  class Color;

  class ColorCreator {
  public:
    ColorCreator(){;}
    ColorCreator(int n){makeColors(n);}
    ~ColorCreator(){purge(colors);}
    
    void addColor(Color* c){colors.push_back(c);}
    int getColorNumber(){return colors.size();};
    Color* getColor(unsigned int i){if (i<colors.size()) return colors[i]; else return 0;}

    void hsvToRgb(double *r, double *g, double *b, double h, double s, double v);

    void makeColors(unsigned int n);

  private:    
    
    /** @link aggregation
     *   @supplierCardinality 1 */
    std::vector <Color*> colors;
  };
}
#endif //COLORCREATOR_H
