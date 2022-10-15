#include "quadtree.hpp"
#include <stdio.h>
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

Quadtree::Quadtree(Rectangle boundary, unsigned int capacity) {
        this->boundary = boundary;
        this->capacity = capacity;
        this->particles = new particle_t[capacity];
        this->particles_count = 0;
        this->northwest = NULL;
        this->northeast = NULL;
        this->southwest = NULL;
        this->southeast = NULL;
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
    this->northwest = new Quadtree(Rectangle(nw_x, nw_y, w/2, h/2), this->capacity);
    this->northeast = new Quadtree(Rectangle(ne_x, ne_y, w/2, h/2), this->capacity);
    this->southwest = new Quadtree(Rectangle(sw_x, sw_y, w/2, h/2), this->capacity);
    this->southeast = new Quadtree(Rectangle(se_x, se_y, w/2, h/2), this->capacity);

    //place the particles in the correct section
    for (int i = 0; i < this->particles_count; i++) {
        if (this->northwest->inboundary(this->particles[i])) {
            this->northwest->insert(this->particles[i]);
            continue;
        }
        else if (this->northeast->inboundary(this->particles[i])) {
            this->northeast->insert(this->particles[i]);
            continue;
        }
        else if (this->southwest->inboundary(this->particles[i])) {
            this->southwest->insert(this->particles[i]);
            continue;
        }
        else if (this->southeast->inboundary(this->particles[i])) {
            this->southeast->insert(this->particles[i]);
        }
    }
    //reasert particles_count to full capacity to prevent future insertion
    this->particles_count = this->capacity;
}

bool Quadtree::inboundary(particle_t particle){
    double x = particle.x;
    double y = particle.y;
    double w = this->boundary.w;
    double h = this->boundary.h;
    double x_min = this->boundary.x - w/2;
    double x_max = this->boundary.x + w/2;
    double y_min = this->boundary.y - h/2;
    double y_max = this->boundary.y + h/2;
    if (x >= x_min && x <= x_max && y >= y_min && y <= y_max) {
        return true;
    }
    return false;
}

void Quadtree::insert(particle_t particle){
    //Do nothing if the particle is not in the boundary
    if (!this->inboundary(particle)) {
        return;
    }

    //if there is space in the quadtree section, add it to the particles array
    if (this->particles_count < this->capacity){
        this->particles[this->particles_count] = particle;
        this->particles_count++;
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
