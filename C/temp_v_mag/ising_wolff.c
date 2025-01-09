#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <time.h>

/*---------------------------------------------------------------------------*\
|                                Parameters                                   |
\*---------------------------------------------------------------------------*/ 

#define N 1000                                   // monte Carlo steps
#define L 250                                   // square lattice size
#define nx L
#define ny L
#define Ti 0.01                                   // temperature
#define Tf 4.0
#define nt 100
#define f 1                                     // saving data frequency

/*---------------------------------------------------------------------------*\
|                              Global variables                               |
\*---------------------------------------------------------------------------*/

int S[nx][ny];          // array of spins
bool C[nx][ny];         // array of clustered spins
int s_add[2][nx * ny];  // array to store cluster positions

// functions
void initial_state();
void random_spin(int *i, int *j);
void cluster_formation(int Si, double p, int i, int j, int *n_add);

int main() {
    int i, j, it, Si, n_add, ic, E, M, ip, jp;
    double T, p;
    FILE *state_file, *mag_file, *energy_file;

    srand(time(NULL));

    // initial state
    initial_state();

    // open files
    mag_file = fopen("data.txt", "w");

    for (int ite = 0; ite <nt; ite++) {

        T = Ti+ite*(Tf-Ti)/nt;
        p = 1.0-exp(-2.0/T);

        M = 0;

        for (it = 1; it <= N; it++) {

            // reset cluster array
            for (i = 0; i < nx; i++) {
                for (j = 0; j < ny; j++) {
                    C[i][j] = false;
                }
            }

            // choose a random spin
            random_spin(&i,&j);
            C[i][j] = true;
            Si = S[i][j];
            n_add = 1;
            s_add[0][0] = i;
            s_add[1][0] = j;

            // form cluster
            for (ic = 0; ic < n_add; ic++) {
                cluster_formation(Si, p, s_add[0][ic], s_add[1][ic], &n_add);
            }

            // flip the cluster
            for (i = 0; i < nx; i++) {
                for (j = 0; j < ny; j++) {
                    if (C[i][j]) {
                        S[i][j] = -S[i][j];
                    }
                }
            }
        }
        M = 0;
        for (i = 0; i < nx; i++) {
            for (j = 0; j < ny; j++) {
            M += S[i][j];
            }
        }
        fprintf(mag_file, "%d %f\n", M, T);
    }

    // close files
    fclose(mag_file);
    return 0;
}

/*---------------------------------------------------------------------------*\
|                                  Functions                                  |
\*---------------------------------------------------------------------------*/

// initial state with random spins
void initial_state() {
    int i, j;
    for (i = 0; i < nx; i++) {
        for (j = 0; j < ny; j++) {
            S[i][j] = (rand()%2)*2-1; // Randomly assign -1 or 1
        }
    }
}

// randomly choose a spin
void random_spin(int *i, int *j) {
    *i = rand()%nx;
    *j = rand()%ny;
}

// form cluster with periodic boundary conditions
void cluster_formation(int Si, double p, int i, int j, int *n_add) {
    int ip, im, jp, jm;
    double r;

    ip = (i+1)%nx;
    im = (i-1+nx)%nx;
    jp = (j+1)%ny;
    jm = (j-1+ny)%ny;

    r = (double)rand()/RAND_MAX;
    if (S[ip][j] == Si && r < p && !C[ip][j]) {
        s_add[0][*n_add] = ip;
        s_add[1][*n_add] = j;
        C[ip][j] = true;
        (*n_add)++;
    }

    r = (double)rand()/RAND_MAX;
    if (S[im][j] == Si && r < p && !C[im][j]) {
        s_add[0][*n_add] = im;
        s_add[1][*n_add] = j;
        C[im][j] = true;
        (*n_add)++;
    }

    r = (double)rand()/RAND_MAX;
    if (S[i][jp] == Si && r < p && !C[i][jp]) {
        s_add[0][*n_add] = i;
        s_add[1][*n_add] = jp;
        C[i][jp] = true;
        (*n_add)++;
    }

    r = (double)rand()/RAND_MAX;
    if (S[i][jm] == Si && r < p && !C[i][jm]) {
        s_add[0][*n_add] = i;
        s_add[1][*n_add] = jm;
        C[i][jm] = true;
        (*n_add)++;
    }
}
