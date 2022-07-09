from tight_binding_approximation import * 
from quspin.operators import hamiltonian 
from quspin.basis import boson_basis_1d 
import numpy as np 

L=L_x*L_y 
alpha=1/5

basis = boson_basis_1d(L,Nb=[0,1,2]) 
# print(basis)

hop=[]
for m in range(L_x*L_y):
    for n in range(L_x*L_y):
        if m in PerBCLat[n]:
            if np.absolute(xy[m][0]-xy[n][0])==L_x-1:
                if xy[m][0] > xy[n][0]:
                    hop.append([-np.exp(-1j*2*np.pi*alpha*xy[m][1]),m,n])
                elif xy[m][0] < xy[n][0]:
                    hop.append([-np.exp(1j*2*np.pi*alpha*xy[m][1]),m,n])
            elif m==n:
                hop.append([0,m,n])
            else:
                if xy[m][0] > xy[n][0]:
                    hop.append([-np.exp(1j*2*np.pi*alpha*xy[m][1]),m,n])
                elif xy[m][0] < xy[n][0]:
                    hop.append([-np.exp(-1j*2*np.pi*alpha*xy[m][1]),m,n])
                else:
                    hop.append([-np.exp(0),m,n])
# print(hop)


U=2

int_1=[]
for m in range(L_x*L_y):
    int_1.append([-U/2,m])
# print(int_1)

#On site interaction
int_2=[]
for m in range(L_x*L_y):
    for n in range(L_x*L_y):
        if m==n:
            int_2.append([U/2,m,n])
# print(int_2)

static=[['+-',hop],['n',int_1],['nn',int_2]]
dynamic=[]
H=hamiltonian(static,dynamic,basis=basis,dtype=np.complex64)
# print(H.todense())

E,V=H.eigh()
