#ifndef GDDCOLORCREATOR_H
#define GDDCOLORCREATOR_H
#include <vector>
#include "detModel/Utilities/GDDpurge.h"

class GDDcolor;

class GDDcolorCreator {
public:
    GDDcolorCreator(){;}
    GDDcolorCreator(int n){makeColors(n);}
    ~GDDcolorCreator(){purge(colors);}

    void addColor(GDDcolor* c){colors.push_back(c);}
    int getColorNumber(){return colors.size();};
    GDDcolor* getColor(int i){if (i<colors.size()) return colors[i]; else return 0;}

    void hsvToRgb(double *r, double *g, double *b, double h, double s, double v);

    void makeColors(int n);

private:    

    /** @link aggregation
         *   @supplierCardinality 1 */
    std::vector <GDDcolor*> colors;
};
#endif //GDDCOLORCREATOR_H
