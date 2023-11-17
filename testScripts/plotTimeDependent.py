import numpy as np
import matplotlib.pyplot as plt
 
 
# Creating dataset
x = np.loadtxt("x.txt")
y = np.loadtxt("y.txt")
Z = np.loadtxt("z.txt")

nx = x.size
ny = y.size

Z = np.resize(Z,[],ny,nx)

[X,Y] = np.meshgrid(x,y)

print(X)
print(Y)
print(Z)

 
# Creating figure
fig = plt.figure(figsize =(14, 9))
ax = plt.axes(projection ='3d')
 
# Creating plot
ax.plot_surface(X, Y, Z[0,:,:])
 
# show plot
plt.show()