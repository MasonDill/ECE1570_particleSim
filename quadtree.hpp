#include "particle.h"
#include <list>
#include <vector> 

class Rectangle {
    public:
        //center point (x,y) and half-width and half-height
        double x, y, w;
        Rectangle();
        Rectangle(double x, double y, double w);
};

class Quadtree {
    public:
        //define the boundary of the quadtree
        Rectangle boundary;
        //define the capacity of a quadtree section 
        unsigned int capacity;
        //models the quadtree section as a single particle at the midpoint of the section
        particle_t* center_of_mass;
        Quadtree* parent;
        //defines the particles in the quadtree section
        //particle_t** particles;

        std::vector<particle_t*> particles; 
        
        void calculateCenterOfMass();
        Quadtree(Rectangle boundary, unsigned int capacity, Quadtree* parent);
        bool insert(particle_t* particle);
        bool hasChildren();
        std::list <Quadtree*>* getLeaves(std::list <Quadtree*>* leaves);
        particle_t* getCenterOfMass();
        bool sharesABorder(Quadtree* other);

    private:
        //defines the four sections of the quadtree
        Quadtree *northwest;
        Quadtree *northeast;
        Quadtree *southwest;
        Quadtree *southeast; 

        void subdivide();
        bool inboundary(particle_t* particle); 
};