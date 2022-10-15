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

Quadtree::Quadtree(Rectangle boundary) {
        this->boundary = boundary;
}

void Quadtree::insert(particle_t particle){
    //TODO
}   
