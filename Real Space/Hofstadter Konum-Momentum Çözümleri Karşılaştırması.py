import numpy as np
# from tight_binding_approximation import *
import matplotlib.pyplot as plt

#Lattice Size
L_x=9
L_y=9

#Create 2-D Lattice. Numbers Represents Each Lattice Sites
lat = np.arange(L_x*L_y).reshape(L_x,L_y)
array = np.append(lat, lat, axis=0)
lattice = array[:-L_x, :]

#Respectively Site Coordinates on Lattice
x_co = np.arange(L_x)
y_co = np.arange(L_y)
arr = []
for j in range(len(x_co)):
    for i in range(len(y_co)):
        arr=np.append(arr, [x_co[i], y_co[j]])
xy = arr.reshape((L_x*L_y,2))

def PerBC(arr):
    neighbors = {}
    for i in range(len(arr)):
        for j, value in enumerate(arr[i]):
            new_neighbors = [
                arr[(i - 1)%L_x][j%L_y],  
                arr[i%L_x][(j + 1)%L_y],  
                arr[(i + 1)%L_x][j%L_y],  
                arr[i%L_x][(j - 1)%L_y]   
            ]
            neighbors[value] = new_neighbors
    return neighbors

PerBCLat = PerBC(lattice)


##################################
# I. KARŞILAŞTIRMA (B=0)
##################################

# def RealH():
#     H = np.zeros((L_x*L_y, L_x*L_y), dtype=complex)
#     for m in range(L_x*L_y):
#         for n in range(L_x*L_y):
#             if m in PerBCLat[n]:
#                 H[m][n] = -1
#     return H

# wR,uR=np.linalg.eig(RealH())
# E_real=np.sort(np.real(wR))

# # kare matris için izinli k değerleri (manyetik alan yok)
# kx=ky=np.array([])
# for m1 in range(L_x):
#     kx=np.append(kx,[-np.pi+(2*np.pi*m1)/L_x])
# for m2 in range(L_y):
#     ky=np.append(ky,[-np.pi+(2*np.pi*m2)/L_y])
# # L_x*L_y tane (kx,ky) çifti var
# E=np.array([])
# for i in kx:
#     for j in ky:
#         E=np.append(E,[2*(np.cos(i)+np.cos(j))],axis=0)
# E_momentum=np.sort(E)

# fig, ax = plt.subplots()
# plt.plot(E_real, 'b-', label=r"$E_\vec{r}$")
# plt.plot(E_momentum, 'r-', label=r"$E_\vec{k}$")
# plt.title("Momentum ve Konum Uzayında Hofstadter Enerjileri")
# ax.set_xlabel(r"$n$")
# ax.set_ylabel(r"$E$")
# plt.legend(prop={'size': 16})
# plt.show()


##################################
# II. KARŞILAŞTIRMA
##################################

q=3 #L_y, q katı
#Hofstadter için izinli k değerleri 
kx=ky=np.array([])
for m1 in range(L_x):
    kx=np.append(kx,[-np.pi+(2*np.pi*m1)/L_x]) #boyutu L_x kadar
for m2 in range(int(L_y/q)):
    ky=np.append(ky,[-np.pi/q+(2*np.pi*m2)/L_y]) #boyutu L_y/q kadar
# L_x*L_y/q tane (kx,ky) çifti var
def MomentumH(alfa, kx, ky):
    M = np.zeros((q,q), dtype=complex)
    for i in range (0, q):
        M[i,i]=2*np.cos(2*np.pi*alfa*i-kx) 
        if i==q-1:
            M[i,i-1]=1
        elif i==0: 
            M[i,i+1]=1
        else: 
            M[i,i-1]=1
            M[i,i+1]=1
        M[0,q-1]=np.exp(q*-1.j*ky)
        M[q-1,0]=np.exp(q*1.j*ky)
    return M
#izinli k değerlerinin enerjileri: (L_x*L_y/q)*q adet öz-değerim var.
E=np.array([])
for i in kx:
    for j in ky:
        wM,uM=np.linalg.eig(MomentumH(1/q, i, j))
        E=np.append(E,wM,axis=0)
E_momentum=np.sort(E)
#izinli k değerleri için bant ayrışımları gözlenir
y = np.zeros(q)
y[:] = 1/q
for i in kx:
    for j in ky:
        x = np.linalg.eigvalsh(MomentumH(1/q, i, j))
        print(x)
        plt.plot(x, y, 'o', markersize=1)

#konum uzayı hofsatdter matrisi        
def RealH(alfa):
    H = np.zeros((L_x*L_y, L_x*L_y), dtype=complex)
    for m in range(L_x*L_y):
        for n in range(L_x*L_y):
            if m in PerBCLat[n]:
                if np.absolute(xy[m][0]-xy[n][0])==L_x-1:
                    if xy[m][0] > xy[n][0]:
                        H[m][n] = -np.exp(-1j*2*np.pi*alfa*xy[m][1])
                    elif xy[m][0] < xy[n][0]:
                        H[m][n] = -np.exp(1j*2*np.pi*alfa*xy[m][1])
                else:
                    if xy[m][0] > xy[n][0]:
                        H[m][n] = -np.exp(1j*2*np.pi*alfa*xy[m][1])
                    elif xy[m][0] < xy[n][0]:
                        H[m][n] = -np.exp(-1j*2*np.pi*alfa*xy[m][1])
                    else:
                        H[m][n] = -np.exp(0)
    return H
#konum uzayı hofstadter matrisi enerjileri
wR,uR=np.linalg.eig(RealH(1/q))
E_real=np.sort(np.real(wR))

fig, ax = plt.subplots()
plt.plot(E_real, 'b-', label=r"$E_\vec{r}$")
plt.plot(E_momentum, 'r-', label=r"$E_\vec{k}$")
plt.title("Momentum ve Konum Uzayında Hofstadter Enerjileri")
ax.set_xlabel(r"$n$")
ax.set_ylabel(r"$E$")
plt.legend(prop={'size': 16})
plt.show()
