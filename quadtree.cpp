#include "quadtree.hpp"
#include <stdio.h>
#define mass 0.01
#define cuttoff 0.01
//
// particle data structure
//

Rectangle::Rectangle() {
    this->x = 0;
    this->y = 0;
    this->w = 0;
}
Rectangle::Rectangle(double x, double y, double w) {
    this->x = x;
    this->y = y;
    this->w = w;
}

Quadtree::Quadtree(Rectangle boundary, unsigned int capacity, Quadtree* parent) {
        this->boundary = boundary;
        this->capacity = capacity;
        // this->particles = new particle_t*[capacity];

        this->northwest = NULL;
        this->northeast = NULL;
        this->southwest = NULL;
        this->southeast = NULL;
        this->parent = parent;
        this->center_of_mass = nullptr;
}

std::list <Quadtree*>* Quadtree::getLeaves(std::list <Quadtree*>* leaves) {
    if(this->hasChildren()) {
        leaves = this->northwest->getLeaves(leaves);
        leaves = this->northeast->getLeaves(leaves);
        leaves = this->southwest->getLeaves(leaves);
        leaves = this->southeast->getLeaves(leaves);
    }
    else {
        if(this->particles.size() > 0)
            leaves->push_back(this);
    }
    return leaves;
}

//divide the quadtree into four sections
void Quadtree::subdivide(){
    double x = this->boundary.x;
    double y = this->boundary.y;
    double w = this->boundary.w;

    double nw_x = x - w/2;
    double nw_y = y + w/2;
    double ne_x = x + w/2;
    double ne_y = y + w/2;
    double sw_x = x - w/2;
    double sw_y = y - w/2;
    double se_x = x + w/2;
    double se_y = y - w/2;
    
    this->northwest = new Quadtree(Rectangle(nw_x, nw_y, w/2), this->capacity, this);
    this->northeast = new Quadtree(Rectangle(ne_x, ne_y, w/2), this->capacity, this);
    this->southwest = new Quadtree(Rectangle(sw_x, sw_y, w/2), this->capacity, this);
    this->southeast = new Quadtree(Rectangle(se_x, se_y, w/2), this->capacity, this);

    //move the particles in the correct subsection, so all particles are always in the leaves
    for (int i = 0; i < this->particles.size(); i++) { //northwest has priority
        this->northwest->insert(this->particles[i]);
        
        this->northeast->insert(this->particles[i]);
        
        this->southwest->insert(this->particles[i]);
        
        this->southeast->insert(this->particles[i]);
            
    }
    //reassert particles_count to full capacity to prevent future insertion
    //this->particles.clear();
}

particle_t* Quadtree::getCenterOfMass() {
    double x = 0;
    double y = 0;
    double m = 0;
    if(this->hasChildren()){
        particle_t nw = *this->northwest->getCenterOfMass();
        particle_t ne = *this->northeast->getCenterOfMass();
        particle_t sw = *this->southwest->getCenterOfMass();
        particle_t se = *this->southeast->getCenterOfMass();

        m = nw.particle_mass + ne.particle_mass + sw.particle_mass + se.particle_mass;
        x = (nw.x * nw.particle_mass + ne.x * ne.particle_mass + sw.x * sw.particle_mass + se.x * se.particle_mass) / (m);
        y = (nw.y * nw.particle_mass + ne.y * ne.particle_mass + sw.y * sw.particle_mass + se.y * se.particle_mass) / (m);
    
        this->center_of_mass = new particle_t;
        this->center_of_mass->x = x;
        this->center_of_mass->y = y;
        this->center_of_mass->particle_mass = m;
    }
    if (this->center_of_mass == nullptr) {
        this->center_of_mass = new particle_t;
        this->center_of_mass->x = 0;
        this->center_of_mass->y = 0;
        this->center_of_mass->particle_mass = 0;
    }
    return this->center_of_mass;
}

bool Quadtree::hasChildren() {
    if (this->northwest == NULL || this->northeast == NULL || this->southwest == NULL || this->southeast == NULL) {
        return false;
    }
    else {
        return true;
    }
}

bool Quadtree::sharesABorder(Quadtree* other) {
    
    if (this->boundary.x == other->boundary.x) {
        if (this->boundary.y - this->boundary.w == other->boundary.y + other->boundary.w || this->boundary.y + this->boundary.w == other->boundary.y - other->boundary.w) {
            return true;
        }
    }
    else if (this->boundary.y == other->boundary.y) {
        if (this->boundary.x - this->boundary.w == other->boundary.x + other->boundary.w || this->boundary.x + this->boundary.w == other->boundary.x - other->boundary.w) {
            return true;
        }
    }
    return false;
}

void Quadtree::calculateCenterOfMass(){
    if(this->center_of_mass != nullptr) {
        delete this->center_of_mass;
    }
    else if (this->particles.size() == 0){
        this->center_of_mass = nullptr;
        return;
    }
    this->center_of_mass = new particle_t;
    this->center_of_mass->x = 0;
    this->center_of_mass->y = 0;
    this->center_of_mass->vx = 0;
    this->center_of_mass->vy = 0;
    this->center_of_mass->particle_mass = mass * this->particles.size();

    for (int i = 0; i < this->particles.size(); i++) {
        //since everything has the same mass, no need to weight the position by mass
        this->center_of_mass->x += this->particles[i]->x;
        this->center_of_mass->y += this->particles[i]->y;
        this->center_of_mass->vx += this->particles[i]->vx;
        this->center_of_mass->vy += this->particles[i]->vy;
    }
    this->center_of_mass->x /= this->particles.size();
    this->center_of_mass->y /= this->particles.size();
    this->center_of_mass->vx /= this->particles.size();
    this->center_of_mass->vy /= this->particles.size();

}

bool Quadtree::inboundary(particle_t* particle){
    double x = particle->x;
    double y = particle->y;

    double w = this->boundary.w;
    double h = this->boundary.w;
    double x_min = this->boundary.x - w - cuttoff;
    double x_max = this->boundary.x + w + cuttoff;
    double y_min = this->boundary.y - h - cuttoff;
    double y_max = this->boundary.y + h + cuttoff;

    //some float math erros can cause the particle to be outside the boundary of all 4 quadrants when its near the border of two quadrants
    if (x >= x_min && x <= x_max && y >= y_min && y <= y_max) {
        return true;
    }
    return false;
}

bool Quadtree::insert(particle_t* particle){
    double C = 1.0;
    double min_width = 0.0001 * C;

    //Do nothing if the particle is not in the boundary
    if (!this->inboundary(particle)) {
        return false;
    }  

    //if there is space in the quadtree section, add it to the particles array
        //TODO CHANGE THIS CONDITIOn
    // if (this->particles.size() < this->capacity || this->boundary.w <= min_width){
    if (this->particles.size() < this->capacity || this->boundary.w <= min_width){
        this->particles.push_back(particle);
        this->calculateCenterOfMass();
        return true;
    }
    //else there is no space in the quadtree section, subdivide it and add the particle to the correct section
    else {
        if (this->northwest == NULL && this->northeast == NULL && this->southwest == NULL && this->southeast == NULL){
            // if(this->boundary.w > min_width)
                this->subdivide();
            // else
            //     this->capacity++;
        }
        
        //add the particle to the correct section
        this->northwest->insert(particle);
        this->northeast->insert(particle);
        this->southwest->insert(particle);
        this->southeast->insert(particle);
    }
    return false;
}   
