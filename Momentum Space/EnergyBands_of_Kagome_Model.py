import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt

pi = np.pi;sin = np.sin;cos = np.cos;sqrt = np.sqrt;exp = np.exp

# t1=1;L1=0;t2=0;L2=0
# t1=1;L1=1;t2=0;L2=0
t1=1;L1=0.28;t2=-0.3;L2=0.2
# t1=1;L1=0.6;t2=0.3;L2=0
def Hamiltonian(k1,k2,k3):
    H = -2*t1*np.array([
    [0, cos(k1), cos(k2)],
    [cos(k1), 0, cos(k3)],        
    [cos(k2), cos(k3), 0]
    ])+2*1j*L1*np.array([
    [0, cos(k1), -cos(k2)],
    [-cos(k1), 0, cos(k3)],        
    [cos(k2), -cos(k3), 0]
    ])-2*t2*np.array([
    [0, cos(k2+k3), cos(k3-k1)],
    [cos(k2+k3), 0, cos(k1+k2)],
    [cos(k3-k1), cos(k1+k2), 0]
    ])+2*1j*L2*np.array([
    [0, -cos(k2+k3), cos(k3-k1)],
    [cos(k2+k3), 0, -cos(k1+k2)],
    [-cos(k3-k1), cos(k1+k2), 0]
    ])
    return H

size = len(Hamiltonian(1, 1, 0))
m=n=25
#infinite lattice
kx_range = np.linspace(np.pi, -np.pi, num=m)
ky_range = np.linspace(np.pi, -np.pi, num=n)
                
energies = np.zeros((m,n,size), dtype=complex) #3D bant için
for i in range(m):
    for j in range(n):
        kx=kx_range[i];ky=ky_range[j]
        k1=kx;k2=kx/2+ky*np.sqrt(3)/2;k3=k2-k1
        H = Hamiltonian(k1, k2, k3)
        evals, evecs = LA.eigh(H)
        energies[i,j,:]=evals
        evalss=np.append(energies, evals) #2D plot için

#3D Plot
X, Y = np.meshgrid(kx_range, ky_range)
fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
# ax.axis([X.min(), X.max(), Y.min(), Y.max()])
for band in range(size):
    ax.plot_surface(X,Y, energies[:,:,band])
ax.set_xlabel(r'$k_{x}$', fontsize=15)
ax.set_ylabel(r'$k_{y}$',fontsize=15)
ax.set_zlabel(r'$E$', fontsize=15)
ax.set_title('Kagome Energy Bands \n'+r'$t_1=$'+str(t1)+','+r'$\lambda_1=$'+str(L1)+'\n'+r'$t_2=$'+str(t2)+','+r'$\lambda_2=$'+str(L2))
ax.view_init(azim=0, elev=0)
plt.show() 
