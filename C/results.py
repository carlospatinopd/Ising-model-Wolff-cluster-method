import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

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
#                          Create the Animation                             #
# ------------------------------------------------------------------------- #

fig, ax = plt.subplots()
im = ax.imshow(Si[:, :, 0], cmap='viridis', vmin=-1, vmax=1)
title = ax.set_title('Spin Configuration at step 0')
ax.set_xlabel('x')
ax.set_ylabel('y')

def update(frame):
    """Update function for the animation."""
    im.set_data(Si[:, :, frame])
    title.set_text(f'Spin Configuration at step {frame * f}')
    return im, title

# Create the animation
ani = FuncAnimation(fig, update, frames=num_snapshots, interval=100, blit=True)

# Save the animation (uncomment one of the following lines based on desired format)
ani.save('state_evolution.mp4', fps=10, dpi=200)

plt.close(fig)  # Close the figure after saving the animation

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
