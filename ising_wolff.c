#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define L 32 // square lattice size
#define T 2.27 // temperature (J/Kb)
#define N 10000 // steps

int S[L][L];

void initial_state() {
    for (int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++) {
            S[i][j] = (rand()%2)*2-1;
        }
    }
}

void pbc(int i) {
    retunr (i+L)%L;
}

int main() {

    double p = 1-exp(-2/T);

    for (int it = 0; it < N; it++) {

    }
    


    srand(time(NULL));

    initial_state();

    double M[N];
    double E[N];

    return 0;
}