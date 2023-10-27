import numpy as np
import matplotlib.pyplot as plt
 
 
# Creating dataset
x = np.loadtxt("x.txt")
y = np.loadtxt("y.txt")
Z = np.loadtxt("z.txt")

nx = x.size
ny = y.size

Z = np.resize(Z,[nx,ny])

[X,Y] = np.meshgrid(y,x)
 
# Creating figure
fig = plt.figure(figsize =(14, 9))
ax = plt.axes(projection ='3d')
 
# Creating plot
ax.plot_surface(X, Y, Z)
 
# show plot
plt.show()