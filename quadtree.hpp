#include "particle.h"
#include <list>
class Rectangle {
    public:
        //center point (x,y) and half-width and half-height
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
        Quadtree* parent;
        //defines the desired width of a quadtree section
        double maximum_interaction_distance;
        unsigned int depth;
        //defines the particles in the quadtree section
        particle_t** particles;

        Quadtree(Rectangle boundary, unsigned int capacity, double max_interaction_distance, unsigned int depth, Quadtree* parent);
        void insert(particle_t* particle);
        bool hasChildren();
        std::list <Quadtree*>* getLeaves(std::list <Quadtree*>* leaves);
        void remove(particle_t* particle);

    private:
        //defines the four sections of the quadtree
        Quadtree *northwest;
        Quadtree *northeast;
        Quadtree *southwest;
        Quadtree *southeast; 

        void subdivide();
        bool inboundary(particle_t* particle); 
};