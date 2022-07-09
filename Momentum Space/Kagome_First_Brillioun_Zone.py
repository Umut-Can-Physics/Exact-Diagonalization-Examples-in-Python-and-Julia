from shapely.geometry import Point, LineString, Polygon, MultiPoint
import random
import numpy as np
import matplotlib.pyplot as plt

# from Kagome_sadece_yakın_komsuluk_matrisi import *
# t1=1;L1 = 0;t2 = 0;L2 = 0:
from Real_Space_Kagome_Model import *

pi = np.pi;sin = np.sin;cos = np.cos;sqrt = np.sqrt;exp = np.exp

N1=N2=5

#Guo ve Franz (Tight-Binding Kagome)
a1=np.array([1,0])
a2=np.array([1/2,sqrt(3)/2])
b1=np.array([-2*pi,2*pi/sqrt(3)])
b2=np.array([0,4*pi/sqrt(3)])

#kagomede dalga vektörü aralıklandırılması
K=np.zeros((0,2))
for m1 in range(-N1,N1):
    for m2 in range(-N2,N2):
        K=np.append(K, [m1/N1*b1+m2/N2*b2], axis=0)  
ind = np.lexsort((K[:,1],K[:,0]))
K = K[ind]

# https://www.researchgate.net/publication/2217014_The_Kagome_Antiferromagnet_A_Schwinger-Boson_Mean-Field_Theory_Study/figures?lo=1
poly = Polygon([(-4*pi/3,0), (-2*pi/3,2*pi/sqrt(3)), (2*pi/3,2*pi/sqrt(3)), (4*pi/3,0), (2*pi/3,-2*pi/sqrt(3)), (-2*pi/3,-2*pi/sqrt(3))])
    
#FBZ içindeki noktalar
k=np.zeros((0,2))
for j in range(len(K)):
    if Point(K[j,0],K[j,1]).within(poly):
        k=np.append(k,np.array([[K[j,0],K[j,1]]]),axis=0)   

#N1=N2=5 iken FBZ içinde nokta seçimi
k = np.delete(k, (0), axis=0)
k = np.delete(k, (3), axis=0)
k=np.append(k,np.array([[K[26,0],K[26,1]]]),axis=0)
k=np.append(k,np.array([[K[14,0],K[14,1]]]),axis=0)
k=np.append(k,np.array([[K[77,0],K[77,1]]]),axis=0)
k=np.append(k,np.array([[K[37,0],K[37,1]]]),axis=0)
#N1=N2=7 iken FBZ içinde nokta seçimi
# k=np.append(k,np.array([[K[34,0],K[34,1]]]),axis=0)
# k=np.append(k,np.array([[K[50,0],K[50,1]]]),axis=0)
# k=np.append(k,np.array([[K[80,0],K[80,1]]]),axis=0)
# k=np.append(k,np.array([[K[109,0],K[109,1]]]),axis=0)
# k=np.append(k,np.array([[K[137,0],K[137,1]]]),axis=0)
# k=np.append(k,np.array([[K[150,0],K[150,1]]]),axis=0)
#N1=N2=4 iken FBZ içinde nokta seçimi (çift sayılarda patlıyoruz)
# k=np.append(k,np.array([[K[12,0],K[12,1]]]),axis=0)
# k=np.append(k,np.array([[K[46,0],K[46,1]]]),axis=0)
# k=np.append(k,np.array([[K[26,0],K[26,1]]]),axis=0)
#N1=N2=3 iken FBZ içinde nokta seçimi
# k=np.append(k,np.array([[K[10,0],K[10,1]]]),axis=0)
# k=np.append(k,np.array([[K[23,0],K[23,1]]]),axis=0)
#N1=N2=6 iken FBZ içinde nokta seçimi
# k=np.append(k,np.array([[K[30,0],K[30,1]]]),axis=0)
# k=np.append(k,np.array([[K[69,0],K[69,1]]]),axis=0)
# k=np.append(k,np.array([[K[105,0],K[105,1]]]),axis=0)
# k=np.append(k,np.array([[K[94,0],K[94,1]]]),axis=0)
# k=np.append(k,np.array([[K[44,0],K[44,1]]]),axis=0)

# tight-binding hamiltonian
def Hamiltonian(k1,k2,k3):
    H = -2*1*np.array([
    [0, np.cos(k1), np.cos(k2)],
    [np.cos(k1), 0, np.cos(k3)],        
    [np.cos(k2), np.cos(k3), 0]
    ])
    return H

#Hangi noktalardan enerjiye nasıl katkılar geliyor?
# kx=K[27,0];ky=K[27,1]
# k1=kx;k2=kx/2+ky*np.sqrt(3)/2;k3=k2-k1
# E,U=np.linalg.eigh(Hamiltonian(k1,k2,k3))
#for eig algorithm and 6x6:
#köşe noktaların enerji katkısı: E=[2,-1,-1] (hangi köşe noktaları aldığımız önemsiz)
#kenar noktaların enerji katkısı: E=[2,-4, 2]

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

#------------------------------------------------------------

#Graphene
# a1=sqrt(3)*np.array([1,0])
# a2=np.array([sqrt(3)/2,3/2])
# b1=2*pi/3*np.array([np.sqrt(3),-1])
# b2=4*pi/3*np.array([0,1])
#Bergholtz ve Liu
# a1=np.array([1,0])
# a2=np.array([1/2,sqrt(3)/2])
# b1=np.array([2*pi,2*pi/sqrt(3)])
# b2=np.array([0,4*pi/sqrt(3)])
#Graphene FBZ
# poly = Polygon([(-4*pi/(3*sqrt(3)),0), (-2*pi/(3*sqrt(3)),2*pi/3), (2*pi/(3*sqrt(3)),2*pi/3), (4*pi/(3*sqrt(3)),0), (2*pi/(3*sqrt(3)),-2*pi/3), (-2*pi/(3*sqrt(3)),-2*pi/3)])
# coord = [[-4*pi/(3*sqrt(3)),0], [-2*pi/(3*sqrt(3)),2*pi/3], [2*pi/(3*sqrt(3)),2*pi/3], [4*pi/(3*sqrt(3)),0], [2*pi/(3*sqrt(3)),-2*pi/3], [-2*pi/(3*sqrt(3)),-2*pi/3]]
# coord = [[-2*pi/3,0], [-pi/3,pi/np.sqrt(3)], [pi/3,pi/np.sqrt(3)], [2*pi/3,0], [pi/3,-pi/np.sqrt(3)], [-pi/3,-pi/np.sqrt(3)]]
#ÇIKTI 1 (iç+köşe+kenar keyfi noktaları) (GRAPHENE için)
#manuel +2 köşe noktası ekleme (N1=N2=3)
# k=np.append(k,np.array([[K[16,0],K[16,1]]]),axis=0)
# k=np.append(k,np.array([[K[8,0],K[8,1]]]),axis=0)
#manuel +2 köşe ve 3 kenar üstünde noktaları ekleme (N1=N2=6)
# k=np.append(k,np.array([[K[56,0],K[56,1]]]),axis=0) #köşe
# k=np.append(k,np.array([[K[106,0],K[106,1]]]),axis=0) #köşe
# k=np.append(k,np.array([[K[81,0],K[81,1]]]),axis=0) 
# k=np.append(k,np.array([[K[42,0],K[42,1]]]),axis=0)
# k=np.append(k,np.array([[K[39,0],K[39,1]]]),axis=0) 
#manuel +2 köşe ve 3 kenar üstünde noktaları ekleme (N1=N2=6)
# k=np.append(k,np.array([[K[50,0],K[50,1]]]),axis=0) #köşe
# k=np.append(k,np.array([[K[100,0],K[100,1]]]),axis=0) #köşe
# k=np.append(k,np.array([[K[117,0],K[117,1]]]),axis=0) 
# k=np.append(k,np.array([[K[114,0],K[114,1]]]),axis=0)
# k=np.append(k,np.array([[K[75,0],K[75,1]]]),axis=0) 
#ÇIKTI 2 (FBZ içine ötelenen nokta)
#FBZ dışında ve içine b1 vektörü ile taşınmış örnek noktalar (N1=N2=7)
# plt.plot(K[91,0],K[91,1],'.', color='blue')
# plt.plot(K[91,0]+b1[0],K[91,1]+b1[1],'.', color='blue')
# plt.plot(K[53,0],K[53,1],'.', color='blue')
# plt.plot(K[53,0]+b1[0],K[53,1]+b1[1],'.', color='blue')
# plt.plot(K[53,0]-b1[0]/6,K[53,1]-b1[1]/6,'.', color='blue')
# plt.plot(K[113,0],K[113,1],'.', color='blue')
# plt.plot(K[126,0]-b1[0],K[126,1]-b1[1],'.', color='blue')
# plt.plot(K[126,0],K[126,1],'.', color='blue')
#(N1=N2=8) 
# plt.plot(K[143,0],K[143,1],'.', color='blue')
# plt.plot(K[143,0]-b2[0],K[143,1]-b2[1],'.', color='blue') 