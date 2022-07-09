from tight_binding_approximation import * 
from quspin.operators import hamiltonian 
from quspin.basis import boson_basis_1d 
import numpy as np 

L=L_x*L_y 
alpha=1/5

basis = boson_basis_1d(L,Nb=[0,1,2]) 
print(basis)

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
print(hop)

static=[['+-',hop]]
dynamic=[]
H=hamiltonian(static,dynamic,basis=basis,dtype=np.complex64)
print(H.todense())

E,V=H.eigh()
