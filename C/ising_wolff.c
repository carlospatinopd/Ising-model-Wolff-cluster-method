#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <time.h>

/*---------------------------------------------------------------------------*\
|                                Parameters                                   |
\*---------------------------------------------------------------------------*/ 

#define N 15000                                  // monte Carlo steps
#define L 200                                   // square lattice size
#define nx L
#define ny L
#define T 3.0                                   // temperature
#define p (1.0-exp(-2.0/T))                     // cluster add probability
#define f 30                                     // saving data frequency

/*---------------------------------------------------------------------------*\
|                              Global variables                               |
\*---------------------------------------------------------------------------*/

int S[nx][ny];          // array of spins
bool C[nx][ny];         // array of clustered spins
int s_add[2][nx * ny];  // array to store cluster positions

// functions
void initial_state();
void random_spin(int *i, int *j);
void cluster_formation(int Si, int i, int j, int *n_add);
void save_data(FILE *state_file, FILE *mag_file, FILE *energy_file, int it);

int main() {
    int i, j, it, Si, n_add, ic, E, M, ip, jp;
    FILE *state_file, *mag_file, *energy_file;

    srand(time(NULL));

    // initial state
    initial_state();

    // open files
    state_file = fopen("state_evolution.txt", "w");
    mag_file = fopen("magnetization.txt", "w");
    energy_file = fopen("energy.txt", "w");

    E = 0;
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
            cluster_formation(Si, s_add[0][ic], s_add[1][ic], &n_add);
        }

        // flip the cluster
        for (i = 0; i < nx; i++) {
            for (j = 0; j < ny; j++) {
                if (C[i][j]) {
                    S[i][j] = -S[i][j];
                }
            }
        }

        // save data
        if (it%f == 0) {
            save_data(state_file,mag_file,energy_file,it);
        }
    }

    // close files
    fclose(state_file);
    fclose(mag_file);
    fclose(energy_file);
    return 0;
}

/*---------------------------------------------------------------------------*\
|                                  Functions                                  |
\*---------------------------------------------------------------------------*/

// initial state with random spins
void initial_state() {
    FILE *param_file, *init_state_file;
    int i, j;

    param_file = fopen("parameters.txt", "w");
    init_state_file = fopen("initial_state.txt", "w");

    fprintf(param_file, "%d %d %d\n", L, N, f);

    for (i = 0; i < nx; i++) {
        for (j = 0; j < ny; j++) {
            S[i][j] = (rand()%2)*2-1; // Randomly assign -1 or 1
            fprintf(init_state_file, "%d\n", S[i][j]);
        }
    }

    fclose(param_file);
    fclose(init_state_file);
}

// randomly choose a spin
void random_spin(int *i, int *j) {
    *i = rand()%nx;
    *j = rand()%ny;
}

// form cluster with periodic boundary conditions
void cluster_formation(int Si, int i, int j, int *n_add) {
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

// save state, magnetization, and energy data each f steps
void save_data(FILE *state_file, FILE *mag_file, FILE *energy_file, int it) {
    int i, j, ip, jp;
    int M = 0, E = 0;

    for (i = 0; i < nx; i++) {
        ip = (i+1)%nx;
        for (j = 0; j < ny; j++) {
            jp = (j+1)%ny;
            M += S[i][j];
            E -= S[i][j]*(S[ip][j]+S[i][jp]);
            fprintf(state_file, "%d\n", S[i][j]);
        }
    }

    fprintf(mag_file, "%d\n", M);
    fprintf(energy_file, "%f\n", E / 2.0);
}
