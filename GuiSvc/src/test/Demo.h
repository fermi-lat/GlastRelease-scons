#ifndef GUI_TEST_DEMO
#define GUI_TEST_DEMO
namespace gui {class DisplayControl; }
class SimpleBox;


class Demo {
public:    
    Demo(gui::DisplayControl& display);
    ~Demo();
    void process();

private:
    int _number;
    SimpleBox* _box1;
    SimpleBox* _box2;
    SimpleBox* _box3;
    SimpleBox* _box4;
    
};



#endif