import numpy as np
import matplotlib.pyplot as plt

# ------------------------------------------------------------------------- #
#                              Load Data                                    #
# ------------------------------------------------------------------------- #

# Load data
parameters = np.loadtxt('parameters.txt')
L = int(parameters[0])  # Lattice size
N = int(parameters[1])  # Total Monte Carlo steps
f = int(parameters[2])  # Frequency of snapshots

magnetization = np.loadtxt('magnetization.txt')
energy = np.loadtxt('energy.txt')
state_evolution = np.loadtxt('state_evolution.txt')

# Initialize 3D array for spin configurations
num_snapshots = N // f
Si = np.zeros((L, L, num_snapshots))

# Populate the 3D spin configuration array
n = 0
for it in range(num_snapshots):
    for i in range(L):
        for j in range(L):
            Si[i, j, it] = state_evolution[n]
            n += 1

# ------------------------------------------------------------------------- #
#                          Plot State Evolution                             #
# ------------------------------------------------------------------------- #

plt.figure(1)
for it in range(num_snapshots):
    plt.imshow(Si[:, :, it], cmap='viridis', vmin=-1, vmax=1)
    plt.colorbar(label='Spin')
    plt.title(f'Spin Configuration at step {it * f}')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.pause(1)
    plt.clf()

# ------------------------------------------------------------------------- #
#                 Plot Magnetization vs Monte Carlo Steps                   #
# ------------------------------------------------------------------------- #

plt.figure(2)
steps = np.arange(1, N + 1, f)
plt.plot(steps, magnetization, linewidth=1.5)
plt.title('Magnetization vs Monte Carlo Steps')
plt.xlabel('Monte Carlo Steps')
plt.ylabel('Magnetization')
plt.grid(True)
plt.show()

# ------------------------------------------------------------------------- #
#                      Plot Energy vs Monte Carlo Steps                     #
# ------------------------------------------------------------------------- #

plt.figure(3)
plt.plot(steps, energy, linewidth=1.5)
plt.title('Energy vs Monte Carlo Steps')
plt.xlabel('Monte Carlo Steps')
plt.ylabel('Energy')
plt.grid(True)
plt.show()
