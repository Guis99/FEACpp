import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d import Axes3D

# Generate some sample data (replace this with your actual data)
# For demonstration purposes, we're creating a sine wave that varies with time.
ns = 50
t = np.linspace(0, 10, ns)
x = np.linspace(-5, 5, ns)
y = np.linspace(-5, 5, ns)
x, y = np.meshgrid(x, y)
data = np.sin(t[:, np.newaxis, np.newaxis] + x**2 + y**2)

# Create a figure and 3D axis
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Set the initial plot for the first frame
plot_surface = ax.plot_surface(x, y, data[0], cmap='viridis')

# Function to update the plot for each frame
def update(frame):
    plot_surface.set_array(data[frame].flatten())
    return plot_surface,

# Create the animation
animation = FuncAnimation(fig, update, frames=len(t), interval=50, blit=True)

# Show the animation
plt.show()
