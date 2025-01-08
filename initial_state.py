import numpy as np
import matplotlib.pyplot as plt

# Parameters
L = 50  # Size of the lattice (from the Fortran code)

# Load data from the file
initial_state = np.loadtxt("initial_state.txt")
final_state = np.loadtxt("final_state.txt")

# Reshape the data into a square lattice
initial_state = initial_state.reshape((L, L))
final_state = final_state.reshape((L, L))

# Plot the lattice
plt.figure(figsize=(16, 8))
plt.subplot(1,2,1)
plt.imshow(initial_state, cmap="tab20b", interpolation="nearest")
plt.title("Initial State of the Ising Model")
plt.xlabel("X-axis")
plt.ylabel("Y-axis")
plt.axis("equal")
plt.subplot(1,2,2)
plt.imshow(final_state, cmap="inferno", interpolation="nearest")
plt.title("Final State of the Ising Model")
plt.xlabel("X-axis")
plt.ylabel("Y-axis")
plt.axis("equal")
plt.show()
