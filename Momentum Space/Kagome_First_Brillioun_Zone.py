from shapely.geometry import Point, LineString, Polygon, MultiPoint
import random
import numpy as np
import matplotlib.pyplot as plt

# t1=1;L1 = 0;t2 = 0;L2 = 0:
from Real_Space_Kagome_Model import *

pi = np.pi;sin = np.sin;cos = np.cos;sqrt = np.sqrt;exp = np.exp

N1=N2=5

a1=np.array([1,0])
a2=np.array([1/2,sqrt(3)/2])
b1=np.array([-2*pi,2*pi/sqrt(3)])
b2=np.array([0,4*pi/sqrt(3)])

# allowed k pairs in kagome
K=np.zeros((0,2))
for m1 in range(-N1,N1):
    for m2 in range(-N2,N2):
        K=np.append(K, [m1/N1*b1+m2/N2*b2], axis=0)  
ind = np.lexsort((K[:,1],K[:,0]))
K = K[ind]

# kagome FBZ
poly = Polygon([(-4*pi/3,0), (-2*pi/3,2*pi/sqrt(3)), (2*pi/3,2*pi/sqrt(3)), (4*pi/3,0), (2*pi/3,-2*pi/sqrt(3)), (-2*pi/3,-2*pi/sqrt(3))])
    
# select points in FBZ
k=np.zeros((0,2))
for j in range(len(K)):
    if Point(K[j,0],K[j,1]).within(poly):
        k=np.append(k,np.array([[K[j,0],K[j,1]]]),axis=0)   

#N1=N2=5
k = np.delete(k, (0), axis=0)
k = np.delete(k, (3), axis=0)
k=np.append(k,np.array([[K[26,0],K[26,1]]]),axis=0)
k=np.append(k,np.array([[K[14,0],K[14,1]]]),axis=0)
k=np.append(k,np.array([[K[77,0],K[77,1]]]),axis=0)
k=np.append(k,np.array([[K[37,0],K[37,1]]]),axis=0)

# tight-binding hamiltonian
def Hamiltonian(k1,k2,k3):
    H = -2*1*np.array([
    [0, np.cos(k1), np.cos(k2)],
    [np.cos(k1), 0, np.cos(k3)],        
    [np.cos(k2), np.cos(k3), 0]
    ])
    return H

E_momentum=np.array([])
for i in range(len(k)):
    kx=k[i,0];ky=k[i,1]
    k1=kx;k2=kx/2+ky*np.sqrt(3)/2;k3=k2-k1
    w,U=np.linalg.eig(Hamiltonian(k1,k2,k3))
    E_momentum=np.append(E_momentum,w,axis=0)
E_momentum=np.real(np.sort(E_momentum))

#PLOT FBZ
fig, axs = plt.subplots()
coord = [[-4*pi/3,0], [-2*pi/3,2*pi/sqrt(3)], [2*pi/3,2*pi/sqrt(3)], [4*pi/3,0], [2*pi/3,-2*pi/sqrt(3)], [-2*pi/3,-2*pi/sqrt(3)]]
order = np.arange(len(K))
for i, txt in enumerate(order):
    axs.annotate(txt, (K[i,0], K[i,1]), fontsize=10)
coord.append(coord[0]) 
xs, ys = zip(*coord) 
axs.plot(xs,ys) 
plt.show() 
axs.plot(K[:,0],K[:,1],'.',color='gray')
axs.scatter(k[:,0],k[:,1], color='red', label="Allowed pairs")  
import matplotlib.patches as patches
plt.annotate(r'$\vec{b}_1$', xy=(0, 0),xytext=(b1[0], b1[1]+0.25),horizontalalignment="center", fontsize=10)
plt.annotate(r'$\vec{b}_2$', xy=(0, 0),xytext=(b2[0], b2[1]),horizontalalignment="center" ,fontsize=10)
style = "Simple, tail_width=1, head_width=5, head_length=5"
kw = dict(arrowstyle=style, color="black")
arrow_origin = np.array([0,0])
arroww=patches.FancyArrowPatch(arrow_origin,np.array([b1[0],b1[1]]), connectionstyle="arc3,rad=0", **kw, label='Reciprocal Vectors')
arrowww=patches.FancyArrowPatch(arrow_origin,np.array([b2[0],b2[1]]), connectionstyle="arc3,rad=0", **kw)
axs.add_patch(arroww)
axs.add_patch(arrowww)
plt.title("Kagome FBZ \n"+str(N1*N2)+r"$\vec{k}$"+" pairs in FBZ")
plt.xlabel(r"$k_x$");plt.ylabel(r"$k_y$")
axs.legend(loc="upper right", scatterpoints=1, fontsize=10)

# PLOT E-n Graph
fig5, ax5 = plt.subplots()
idx = np.argsort(E_momentum)
sorted = E_momentum[idx]
x_sort = np.arange(0, len(E_momentum), 1)
y_sort = sorted
ax5.set_title("Energies of Kagome Model in Momentum and Real Spaces \n"+str(N1)+r"$\times $"+str(N2)+r"$\quad (t_1=1,t_2=\lambda_1=\lambda_2=0)$")
ax5.set_xlabel(r"$n$")
ax5.set_ylabel(r"$E$")
ax5.plot(x_sort, y_sort, 'b-', markersize=0.7, label=r"$E_\vec{k}$")
ax5.plot(x_sort_real, y_sort_real, 'r.', markersize=5, label=r"$E_\vec{r}$")
plt.legend(prop={'size': 16}, scatterpoints=15)
