import numpy as np
import matplotlib.pyplot as plt

# ------------------------------------------------------------------------- #
#                              Load Data                                    #
# ------------------------------------------------------------------------- #

parameters = np.loadtxt('parameters.txt')
L = int(parameters[0])
N = int(parameters[1])
f = int(parameters[2])
magnetization = np.loadtxt('magnetization.txt')
energy = np.loadtxt('energy.txt')
state_evolution = np.loadtxt('state_evolution.txt')
num_snapshots = N//f
Si = np.zeros((L, L, num_snapshots))
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
    plt.title(f'Spin Configuration at step {it*f}')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.pause(0.01)
    plt.clf()

# ------------------------------------------------------------------------- #
#                 Plot Magnetization vs Monte Carlo Steps                   #
# ------------------------------------------------------------------------- #

plt.figure(2)
steps = np.arange(1, N + 1, f)
plt.plot(steps, magnetization, linewidth=1.5)
plt.xlabel('Monte Carlo Steps')
plt.ylabel('Magnetization')
plt.show()

# ------------------------------------------------------------------------- #
#                      Plot Energy vs Monte Carlo Steps                     #
# ------------------------------------------------------------------------- #

plt.figure(3)
plt.plot(steps, energy, linewidth=1.5)
plt.xlabel('Monte Carlo Steps')
plt.ylabel('Energy')
plt.show()
