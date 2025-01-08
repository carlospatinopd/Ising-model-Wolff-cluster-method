#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <random>

const int N = 100000;  // Monte Carlo steps
const int L = 50;      // Square lattice size
const double T = 2.2;  // Temperature
const double p = 1.0 - exp(-2.0 / T);  // Bond probability

std::mt19937 rng(std::random_device{}());
std::uniform_real_distribution<> random_double(0.0, 1.0);

using Matrix = std::vector<std::vector<int>>;
using BoolMatrix = std::vector<std::vector<bool>>;

// Function declarations
void initial_state(Matrix &S);
void random_spin(int &i, int &j);
void cluster_formation(int i, int j, int Si, int &n_add, std::vector<std::pair<int, int>> &s_add, BoolMatrix &C, const Matrix &S);
void save_to_file(const Matrix &S, const std::string &filename);

int main() {
    Matrix S(L, std::vector<int>(L));  // Array of spins
    BoolMatrix C(L, std::vector<bool>(L, false));  // Clustered spins
    std::vector<std::pair<int, int>> s_add(4);  // Neighbors for cluster formation

    initial_state(S);
    save_to_file(S, "initial_state.txt");

    for (int it = 0; it < N; ++it) {
        std::fill(C.begin(), C.end(), std::vector<bool>(L, false));

        int i, j, Si;
        random_spin(i, j);
        C[i][j] = true;
        Si = S[i][j];

        int n_add = 0;
        cluster_formation(i, j, Si, n_add, s_add, C, S);

        while (n_add > 0) {
            std::vector<std::pair<int, int>> new_s_add;
            for (int ic = 0; ic < n_add; ++ic) {
                int x = s_add[ic].first;
                int y = s_add[ic].second;
                cluster_formation(x, y, Si, n_add, new_s_add, C, S);
            }
            s_add = new_s_add;
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
    return 0;
}

// Initialize spins randomly (-1 or 1)
void initial_state(Matrix &S) {
    for (int i = 0; i < L; ++i) {
        for (int j = 0; j < L; ++j) {
            S[i][j] = (random_double(rng) < 0.5) ? -1 : 1;
        }
    }
}

// Choose a random spin
void random_spin(int &i, int &j) {
    i = rng() % L;
    j = rng() % L;
}

// Cluster formation with periodic boundary conditions
void cluster_formation(int i, int j, int Si, int &n_add, std::vector<std::pair<int, int>> &s_add, BoolMatrix &C, const Matrix &S) {
    int ip = (i + 1) % L, im = (i - 1 + L) % L;
    int jp = (j + 1) % L, jm = (j - 1 + L) % L;

    n_add = 0;

    if (S[ip][j] == Si && random_double(rng) < p && !C[ip][j]) {
        s_add.emplace_back(ip, j);
        C[ip][j] = true;
        ++n_add;
    }
    if (S[im][j] == Si && random_double(rng) < p && !C[im][j]) {
        s_add.emplace_back(im, j);
        C[im][j] = true;
        ++n_add;
    }
    if (S[i][jp] == Si && random_double(rng) < p && !C[i][jp]) {
        s_add.emplace_back(i, jp);
        C[i][jp] = true;
        ++n_add;
    }
    if (S[i][jm] == Si && random_double(rng) < p && !C[i][jm]) {
        s_add.emplace_back(i, jm);
        C[i][jm] = true;
        ++n_add;
    }
}

// Save lattice state to file
void save_to_file(const Matrix &S, const std::string &filename) {
    std::ofstream file(filename);
    for (const auto &row : S) {
        for (const auto &val : row) {
            file << val << " ";
        }
        file << "\n";
    }
}
