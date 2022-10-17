#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <list>
#include "common.h"
#include "quadtree.hpp"
//
//  benchmarking program

int main( int argc, char **argv )
{    
    int navg,nabsavg=0;
    double davg,dmin, absmin=1.0, absavg=0.0;

    if( find_option( argc, argv, "-h" ) >= 0 )
    {
        printf( "Options:\n" );
        printf( "-h to see this help\n" );
        printf( "-n <int> to set the number of particles\n" );
        printf( "-o <filename> to specify the output file name\n" );
        printf( "-s <filename> to specify a summary file name\n" );
        printf( "-c <int> to set the capcity of the quadtree\n" );
        printf( "-no turns off all correctness checks and particle output\n");
        return 0;
    }
    //16 is the default capacity of the quadtree
    int n = read_int( argc, argv, "-n", 1000 );
    int capacity = read_int( argc, argv, "-c", 100);

    char *savename = read_string( argc, argv, "-o", NULL );
    char *sumname = read_string( argc, argv, "-s", NULL );
    
    FILE *fsave = savename ? fopen( savename, "w" ) : NULL;
    FILE *fsum = sumname ? fopen ( sumname, "a" ) : NULL;

    //
    //initialize the particles
    //
    particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );
    set_size( n );
    init_particles( n, particles );

    //
    // initialize the quadtree
    //

    double middle = get_size() * 0.5;
    Rectangle boundary = Rectangle(middle, middle, middle);

    //
    //  simulate a number of time steps
    //
    double simulation_time = read_timer( );
    
    for( int step = 0; step < NSTEPS; step++ )
    {
	navg = 0;
    davg = 0.0;
	dmin = 1.0; //dmin = 1.0 -> r2 always > cuttoff^2, no particles are interacting 

        //create the quadtree every time step to account for movement of particles
        Quadtree* tree = new Quadtree(boundary, capacity, nullptr);
        for (int i = 0; i < n; i++) {
            if (!tree->insert(&particles[i]))
                printf("Error: particle could not be inserted into the quadtree");
        }

        for( int i = 0; i < n; i++ )
        {
            particles[i].ax = particles[i].ay = 0;
        }

        std::list <Quadtree*>* leaves = tree->getLeaves(new std::list <Quadtree*>());
        for (std::list<Quadtree*>::iterator it = leaves->begin(); it != leaves->end(); ++it) {
            Quadtree* subquad = *it;
            //  for each particle in the subquadtree section,
            for (int i = 0; i < subquad->particles.size(); i++) {
                //calculate the forces of the particles on each other in the subquadtree
                for (int j = 0; j < subquad->particles.size(); j++) {
                    if(i != j)
                        apply_force( *subquad->particles[i], *subquad->particles[j],&dmin,&davg,&navg);
                }
            }
            for(std::list<Quadtree*>::iterator it2 = leaves->begin(); it2 != leaves->end(); ++it2) {
                Quadtree* quad2 = *it2;
                //  for each particle in the subquadtree section,
                for (int i = 0; i < subquad->particles.size(); i++) {
                    //calculate the forces of the particles on each other in the subquadtree
                    for (int j = 0; j < quad2->particles.size(); j++) {
                        apply_force( *subquad->particles[i], *quad2->particles[j],&dmin,&davg,&navg);
                    }
                }
            }
        }


        //  move particles
        for( int i = 0; i < n; i++ ) 
            move( particles[i] );		

        if( find_option( argc, argv, "-no" ) == -1 )
        {
          // Computing statistical data
          if (navg) {
            absavg +=  davg/navg;
            nabsavg++;
          }

          if (dmin < absmin)
            absmin = dmin;
    
          //  save if necessary
          if( fsave && (step%SAVEFREQ) == 0 )
              save( fsave, n, particles );
        }
    }
    simulation_time = read_timer( ) - simulation_time;
    
    printf( "n = %d, simulation time = %g seconds", n, simulation_time);

    if( find_option( argc, argv, "-no" ) == -1 )
    {
      if (nabsavg) absavg /= nabsavg;
    // 
    //  -The minimum distance absmin between 2 particles during the run of the simulation
    //  -A Correct simulation will have particles stay at greater than 0.4 (of cutoff) with typical values between .7-.8
    //  -A simulation where particles don't interact correctly will be less than 0.4 (of cutoff) with typical values between .01-.05
    //
    //  -The average distance absavg is ~.95 when most particles are interacting correctly and ~.66 when no particles are interacting
    //
    printf( ", absmin = %lf, absavg = %lf", absmin, absavg);
    if (absmin < 0.4) printf ("\nThe minimum distance is below 0.4 meaning that some particle is not interacting");
    if (absavg < 0.8) printf ("\nThe average distance is below 0.8 meaning that most particles are not interacting");
    }
    printf("\n");

    //
    // Printing summary data
    //
    if( fsum) 
        fprintf(fsum,"%d %g\n",n,simulation_time);
 
    //
    // Clearing space
    //
    if( fsum )
        fclose( fsum );    
    free( particles );
    if( fsave )
        fclose( fsave );
    
    return 0;
}
