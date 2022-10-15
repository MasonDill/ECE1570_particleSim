#include "particle.h"
#include <list>
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
        //defines the number of particles in the quadtree section
        unsigned int particles_count;
        //models the quadtree section as a single particle at the midpoint of the section
        particle_t center_of_mass;
        //defines the particles in the quadtree section
        particle_t** particles;

        Quadtree(Rectangle boundary, unsigned int capacity);
        void insert(particle_t* particle);
        bool hasChildren();
        std::list <Quadtree*>* getLeaves(std::list <Quadtree*>* leaves);

    private:
        //defines the four sections of the quadtree
        Quadtree *northwest;
        Quadtree *northeast;
        Quadtree *southwest;
        Quadtree *southeast; 

        void subdivide();
        bool inboundary(particle_t* particle); 
};