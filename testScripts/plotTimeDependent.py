import numpy as np
import matplotlib.pyplot as plt
 
 
# Creating dataset
x = np.loadtxt("xt.txt")
y = np.loadtxt("yt.txt")
Z = np.loadtxt("zt.txt")
print(Z.shape)

nx = x.size
ny = y.size

Z = Z.reshape(-1,ny,nx)

print(Z.shape)

[X,Y] = np.meshgrid(x,y)

# print(X)
# print(Y)
# print(Z)

 
# Creating figure
fig = plt.figure(figsize =(14, 9))
ax = plt.axes(projection ='3d')
ax.set_zlim([0,1])
# Creating plot
ax.plot_surface(X, Y, Z[4])

 
# show plot
plt.show()