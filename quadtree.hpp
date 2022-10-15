#include "particle.h"
class Rectangle {
    double x, y, w, h;
    public:
        Rectangle();
        Rectangle(double x, double y, double w, double h);
};

class Quadtree {
    Rectangle boundary;
    public:
        Quadtree(Rectangle boundary);
        void insert(particle_t particle);     
};