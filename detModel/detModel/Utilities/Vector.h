#ifndef VECTOR_H
#define VECTOR_H
#define GDDPI 3.1415926536

namespace detModel{

  class Vector{
  public:
    
    Vector():x(0),y(0),z(0){};
    Vector(double px, double py, double pz):x(px), y(py), z(pz){};

    double getX(){return x;}
    double getY(){return y;}
    double getZ(){return z;}

    void setX(float a){x = a;}
    void setY(float a){y = a;}
    void setZ(float a){z = a;}

    void Equate(Vector* v);
    void Add(Vector* t);
    void Multiply(float a);
    void Rotate(int dir, double theta);
  

    double x;
    double y;
    double z;

  };
}
#endif // VECTORH


