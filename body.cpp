#include <iostream>
#include <string>
#include <vector>
#include <math.h>
#include <stdio.h>

#define N 3 // # of particles
#define s 1e-5 // softening
#define G 4.4988e-3 // pc^3 Msun^-1 Myr^-2

// custom data type
typedef struct {
    double mass;
    double pos[3]; // x,y,z coords
    double vel[3];
    double acc[3];
} Particle;

void compute_forces(Particle *p){
    for (int i = 0; i < N; ++i){
        // zero out acceleration vector for particle i
        p[i].acc[0] = p[i].acc[1] = p[i].acc[2] = 0.0;

        for (int j = 0; j < N; ++j){
            if (i ==j) continue; // skip gravitational force on itself

            // displacement vector
            double dx = p[j].pos[0] - p[i].pos[0];
            double dy = p[j].pos[1] - p[i].pos[1];
            double dz = p[j].pos[2] - p[i].pos[2];
            
            // softened squared and cubed distances
            double dist_sq = dx*dx + dy*dy + dz*dz + s*s;
            double dist_3 = pow(dist_sq, 1.5);

            // scalar acceleration due to grav
            double force = G * p[j].mass / dist_3;

            // project scalar onto vector components
            p[i].acc[0] += dx * force;
            p[i].acc[1] += dy * force;
            p[i].acc[2] += dz * force;

        }
    }

}

int main() {

    // initialize particles
    Particle p[N] = {
        {1.0,  {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}},
        {2.0,  {1.0, 0.0, 0.0}, {0.0, 0.0, 0.0}},
        {1.5,  {0.0, 1.0, 0.0}, {0.0, 0.0, 0.0}}
    };

    compute_forces(p);

    for (int i = 0; i < N; ++i){
        printf("Particle %d: acc = (%.5e, %.5e, %.5e)\n", 
            i, p[i].acc[0], p[i].acc[1], p[i].acc[2]);
    }

    return 0;
}