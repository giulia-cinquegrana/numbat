#include <iostream>
#include <math.h>
#include <stdio.h>

#define s 1e-5 // softening
#define G 4.4988e-3 // pc^3 Msun^-1 Myr^-2

// custom data type
typedef struct {
    double mass;
    double pos[3]; // x,y,z coords
    double vel[3];
    double acc[3];
} Particle;

void initialize_particles(Particle* p, int N) {

    srand(time(NULL)); // each run gets different ICs

    for (int i = 0; i < N; ++i){
        // random mass between 0.5 and 2.0 Msun
        p[i].mass = 0.5 + (rand() / (double)RAND_MAX) * 1.5;
        // disk geometry
        // random angle from 0 to 2pi
        double angle = 2.0 * M_PI * (rand() / (double)RAND_MAX);
        // random radius, 
        // weighted to form higher densities toward centre
        double radius = 5.0 * sqrt(rand() / (double)RAND_MAX);
        // particles in flat disk
        p[i].pos[0] = radius * cos(angle);
        p[i].pos[1] = radius * sin(angle);
        p[i].pos[2] = 0.0;

        // circular motion
        double v = sqrt(G * p[i].mass / (radius + s));

        p[i].vel[0] = -v * sin(angle);
        p[i].vel[1] =  v * cos(angle);
        p[i].vel[2] = 0.0;

    }


}

void compute_forces(Particle *p, int N){

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

    int N = 100;

    double t = 0.0;
    double dt = 0.01;
    int nsteps = 1000;

    Particle* p = new Particle[N];

    initialize_particles(p, N);

    compute_forces(p, N);

    // leapfrog. kick, drift, re- compute_forces, kick

    for (int step=0; step < nsteps; ++step){
        for (int i = 0; i < N; ++i){
            for (int j = 0; j < 3; ++j){
                p[i].vel[j]+= 0.5 * dt * p[i].acc[j];
            }

        }
        for (int i = 0; i < N; ++i){
            for (int j = 0; j < 3; ++j){
                p[i].pos[j] += dt * p[i].vel[j];
            }
        }
        compute_forces(p, N);
        
        for (int i = 0; i < N; ++i){
            for (int j = 0; j < 3; ++j){
                p[i].vel[j]+= 0.5 * dt * p[i].acc[j];
            }
        }
    }

    // for (int i = 0; i < N; ++i){
    //     printf("Particle %d: acc = (%.5e, %.5e, %.5e)\n", 
    //         i, p[i].acc[0], p[i].acc[1], p[i].acc[2]);
    // }


    delete[] p;
    return 0;
}