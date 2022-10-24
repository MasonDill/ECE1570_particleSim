#ifndef PARTICLE_H
#define PARTICLE_H
typedef struct 
{
  double x;
  double y;
  double vx;
  double vy;
  double ax;
  double ay;
  double particle_mass;
} particle_t;
#endif