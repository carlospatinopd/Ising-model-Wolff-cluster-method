#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <time.h>

#define N 100000  // Monte Carlo steps
#define L 50      // Square lattice size
#define T 2.2     // Temperature
#define P (1.0 - exp(-2.0 / T))  // Bond probability

// Function declarations
void initial_state(int **S);
void random_spin(int *i, int *j);
void cluster_formation(int i, int j, int Si, int *n_add, int s_add[][2], bool **C, int **S);
void save_to_file(int **S, const char *filename);
int periodic_boundary(int index, int limit);

int main() {
    // Allocate 2D arrays
    int **S = malloc(L * sizeof(int *));
    bool **C = malloc(L * sizeof(bool *));
    for (int i = 0; i < L; i++) {
        S[i] = malloc(L * sizeof(int));
        C[i] = malloc(L * sizeof(bool));
    }

    int s_add[4][2];  // Neighbors for cluster formation
    srand((unsigned int)time(NULL));  // Seed random number generator

    initial_state(S);
    save_to_file(S, "initial_state.txt");

    for (int it = 0; it < N; ++it) {
        for (int i = 0; i < L; ++i)
            for (int j = 0; j < L; ++j)
                C[i][j] = false;

        int i, j, Si;
        random_spin(&i, &j);
        C[i][j] = true;
        Si = S[i][j];

        int n_add = 0;
        cluster_formation(i, j, Si, &n_add, s_add, C, S);

        while (n_add > 0) {
            int new_s_add[4][2];
            int new_n_add = 0;

            for (int ic = 0; ic < n_add; ++ic) {
                int x = s_add[ic][0];
                int y = s_add[ic][1];
                cluster_formation(x, y, Si, &new_n_add, new_s_add, C, S);
            }

            for (int k = 0; k < new_n_add; ++k) {
                s_add[k][0] = new_s_add[k][0];
                s_add[k][1] = new_s_add[k][1];
            }
            n_add = new_n_add;
        }

        for (int x = 0; x < L; ++x) {
            for (int y = 0; y < L; ++y) {
                if (C[x][y]) {
                    S[x][y] = -S[x][y];
                }
            }
        }
    }

    save_to_file(S, "final_state.txt");

    // Free allocated memory
    for (int i = 0; i < L; i++) {
        free(S[i]);
        free(C[i]);
    }
    free(S);
    free(C);

    return 0;
}

// Initialize spins randomly (-1 or 1)
void initial_state(int **S) {
    for (int i = 0; i < L; ++i) {
        for (int j = 0; j < L; ++j) {
            S[i][j] = (rand() / (double)RAND_MAX < 0.5) ? -1 : 1;
        }
    }
}

// Choose a random spin
void random_spin(int *i, int *j) {
    *i = rand() % L;
    *j = rand() % L;
}

// Cluster formation with periodic boundary conditions
void cluster_formation(int i, int j, int Si, int *n_add, int s_add[][2], bool **C, int **S) {
    int ip = periodic_boundary(i + 1, L);
    int im = periodic_boundary(i - 1, L);
    int jp = periodic_boundary(j + 1, L);
    int jm = periodic_boundary(j - 1, L);

    *n_add = 0;

    if (S[ip][j] == Si && (rand() / (double)RAND_MAX) < P && !C[ip][j]) {
        s_add[*n_add][0] = ip;
        s_add[*n_add][1] = j;
        C[ip][j] = true;
        (*n_add)++;
    }
    if (S[im][j] == Si && (rand() / (double)RAND_MAX) < P && !C[im][j]) {
        s_add[*n_add][0] = im;
        s_add[*n_add][1] = j;
        C[im][j] = true;
        (*n_add)++;
    }
    if (S[i][jp] == Si && (rand() / (double)RAND_MAX) < P && !C[i][jp]) {
        s_add[*n_add][0] = i;
        s_add[*n_add][1] = jp;
        C[i][jp] = true;
        (*n_add)++;
    }
    if (S[i][jm] == Si && (rand() / (double)RAND_MAX) < P && !C[i][jm]) {
        s_add[*n_add][0] = i;
        s_add[*n_add][1] = jm;
        C[i][jm] = true;
        (*n_add)++;
    }
}

// Save lattice state to file
void save_to_file(int **S, const char *filename) {
    FILE *file = fopen(filename, "w");
    if (!file) {
        perror("Error opening file");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < L; ++i) {
        for (int j = 0; j < L; ++j) {
            fprintf(file, "%d ", S[i][j]);
        }
        fprintf(file, "\n");
    }

    fclose(file);
}

// Apply periodic boundary conditions
int periodic_boundary(int index, int limit) {
    if (index < 0) return limit - 1;
    if (index >= limit) return 0;
    return index;
}
