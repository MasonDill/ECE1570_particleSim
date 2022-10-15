#include "particle.h"
class Rectangle {
    public:
        double x, y, w, h;
        Rectangle();
        Rectangle(double x, double y, double w, double h);
};

class Quadtree {
    public:
        //define the boundary of the quadtree
        Rectangle boundary;
        //define the capacity of a quadtree section 
        unsigned int capacity;
        //defines the particles in the quadtree section
        particle_t *particles;
        //models the quadtree section as a single particle
        particle_t center;

        Quadtree(Rectangle boundary, unsigned int capacity);
        void insert(particle_t particle);
    private:
        //defines the number of particles in the quadtree section
        unsigned int particles_count;

        Quadtree *northwest;
        Quadtree *northeast;
        Quadtree *southwest;
        Quadtree *southeast; 

        void subdivide();
        bool inboundary(particle_t particle);
           
};