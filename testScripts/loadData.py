import numpy as np

X = np.loadtxt("x.txt")
Y = np.loadtxt("y.txt")
Z = np.loadtxt("z.txt")

nx = X.size
ny = Y.size

Z = np.resize(Z,[nx,ny])

nz = Z.shape

print(X)
print(Y)
print(Z)
print(nx,ny,nz)