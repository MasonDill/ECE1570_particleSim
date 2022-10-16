#include "quadtree.hpp"
#include <stdio.h>
#define mass 0.01
//
// particle data structure
//

Rectangle::Rectangle() {
    this->x = 0;
    this->y = 0;
    this->w = 0;
    this->h = 0;
}
Rectangle::Rectangle(double x, double y, double w, double h) {
    this->x = x;
    this->y = y;
    this->w = w;
    this->h = h;
}

Quadtree::Quadtree(Rectangle boundary, unsigned int capacity, Quadtree* parent) {
        this->boundary = boundary;
        this->capacity = capacity;
        this->particles = new particle_t*[capacity];
        this->particles_count = 0;

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
        leaves->push_back(this);
    }
    return leaves;
}

//divide the quadtree into four sections
void Quadtree::subdivide(){
    double x = this->boundary.x;
    double y = this->boundary.y;
    double w = this->boundary.w;
    double h = this->boundary.h;

    double nw_x = x - w/2;
    double nw_y = y + h/2;
    double ne_x = x + w/2;
    double ne_y = y + h/2;
    double sw_x = x - w/2;
    double sw_y = y - h/2;
    double se_x = x + w/2;
    double se_y = y - h/2;
    
    this->northwest = new Quadtree(Rectangle(nw_x, nw_y, w/2, h/2), this->capacity, this);
    this->northeast = new Quadtree(Rectangle(ne_x, ne_y, w/2, h/2), this->capacity, this);
    this->southwest = new Quadtree(Rectangle(sw_x, sw_y, w/2, h/2), this->capacity, this);
    this->southeast = new Quadtree(Rectangle(se_x, se_y, w/2, h/2), this->capacity, this);

    //place the particles in the correct section
    for (int i = 0; i < this->particles_count; i++) {
        if (this->northwest->inboundary(this->particles[i])) {
            this->northwest->insert(this->particles[i]);
            this->particles[i] = NULL;
        }
        else if (this->northeast->inboundary(this->particles[i])) {
            this->northeast->insert(this->particles[i]);
            this->particles[i] = NULL;
        }
        else if (this->southwest->inboundary(this->particles[i])) {
            this->southwest->insert(this->particles[i]);
            this->particles[i] = NULL;
        }
        else if (this->southeast->inboundary(this->particles[i])) {
            this->southeast->insert(this->particles[i]);
            this->particles[i] = NULL;
        }
        else
            printf("Error: particle not in any section");
    }
    //reassert particles_count to full capacity to prevent future insertion
    this->particles_count = this->capacity;
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
    if (this->boundary.x == other->boundary.x && this->boundary.y == other->boundary.y) {
        return false;
    }
    else if (this->boundary.x == other->boundary.x) {
        if (this->boundary.y + this->boundary.h/2 == other->boundary.y - other->boundary.h/2) {
            return true;
        }
        else if (this->boundary.y - this->boundary.h/2 == other->boundary.y + other->boundary.h/2) {
            return true;
        }
        else {
            return false;
        }
    }
    else if (this->boundary.y == other->boundary.y) {
        if (this->boundary.x + this->boundary.w/2 == other->boundary.x - other->boundary.w/2) {
            return true;
        }
        else if (this->boundary.x - this->boundary.w/2 == other->boundary.x + other->boundary.w/2) {
            return true;
        }
        else {
            return false;
        }
    }
    else {
        return false;
    }
}

void Quadtree::calculateCenterOfMass(){
    if(this->center_of_mass != nullptr) {
        delete this->center_of_mass;
    }
    this->center_of_mass = new particle_t;
    this->center_of_mass->x = 0;
    this->center_of_mass->y = 0;
    this->center_of_mass->vx = 0;
    this->center_of_mass->vy = 0;
    this->center_of_mass->particle_mass = mass * this->particles_count;

    for (int i = 0; i < this->particles_count; i++) {
        //since everything has the same mass, no need to weight the position by mass
        this->center_of_mass->x += this->particles[i]->x;
        this->center_of_mass->y += this->particles[i]->y;
        this->center_of_mass->vx += this->particles[i]->vx;
        this->center_of_mass->vy += this->particles[i]->vy;
    }
    this->center_of_mass->x /= this->particles_count;
    this->center_of_mass->y /= this->particles_count;
    this->center_of_mass->vx /= this->particles_count;
    this->center_of_mass->vy /= this->particles_count;

}

bool Quadtree::inboundary(particle_t* particle){
    double x = particle->x;
    double y = particle->y;

    double w = this->boundary.w;
    double h = this->boundary.h;
    double x_min = this->boundary.x - w;
    double x_max = this->boundary.x + w;
    double y_min = this->boundary.y - h;
    double y_max = this->boundary.y + h;
    if (x >= x_min && x <= x_max && y >= y_min && y <= y_max) {
        return true;
    }
    return false;
}

void Quadtree::insert(particle_t* particle){
    //Do nothing if the particle is not in the boundary
    if (!this->inboundary(particle)) {
        return;
    }

    //if there is space in the quadtree section, add it to the particles array
    if (this->particles_count < this->capacity){
        this->particles[this->particles_count] = particle;
        this->particles_count++;
        this->calculateCenterOfMass();
        return;
    }
    //else there is no space in the quadtree section, subdivide it and add the particle to the correct section
    else {
        if (this->northwest == NULL && this->northeast == NULL && this->southwest == NULL && this->southeast == NULL){
            this->subdivide();
        }
        //add the particle to the correct section
        this->northwest->insert(particle);
        this->northeast->insert(particle);
        this->southwest->insert(particle);
        this->southeast->insert(particle);
    }

}   
