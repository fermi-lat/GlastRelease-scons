#ifndef COLOR_H
#define COLOR_H

namespace detModel{

  class Color {
  public:
    
    Color():red(0),green(0),blue(0),tra(0){;}
    Color(float r, float g, float b, float t):red(r),green(g),blue(b),tra(t){;}
    
    float getTra(){return tra;};
    void setTra(float t){tra = t;};
    float getBlue(){return blue;};
    void setBlue(float b){blue = b;};
    float getGreen(){return green;};
    void setGreen(float g){green = g;};
    float getRed(){return red;};
    void setRed(float r){red = r;};
    
  private:
    
    float red;
    float green;
    float blue;
    float tra;
  };
}
#endif //COLOR_H
