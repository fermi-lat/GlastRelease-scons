#ifndef GDDMATRIX_H
#define GDDMATRIX_H

#define GDDPI 3.1415926536

class GDDvector{
 public:

  GDDvector():x(0),y(0),z(0){};
  GDDvector(double px, double py, double pz):x(px), y(py), z(pz){};

  double getX(){return x;}
  double getY(){return y;}
  double getZ(){return z;}

  void setX(float a){x = a;}
  void setY(float a){y = a;}
  void setZ(float a){z = a;}

  double x;
  double y;
  double z;
};

void GDDequate(GDDvector* t, GDDvector* v);
void GDDtranslate(GDDvector* t, GDDvector* v);
void GDDrotate(int dir, double theta, GDDvector* v);


#endif // GDDMATRIX_H
