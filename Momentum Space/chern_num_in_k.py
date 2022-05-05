import numpy as np

size=N=5
t = 1
p = 1
q = size 
phi = p/q

def HMat(phi, kx, ky):
    Hk = np.zeros((size,size), dtype=complex)
    for i in range(0, size):
        Hk[i,i] = -2*t*np.cos(ky-2*i*np.pi*phi)
        if i==size-1:
            Hk[i,i-1] = 1
        elif i==0:
            Hk[i,i+1] = 1
        else:
            Hk[i,i-1] = 1
            Hk[i,i+1] = 1
        Hk[0,size-1]= -t*np.exp(-q*1.j*kx)
        Hk[size-1,0]= -t*np.exp(q*1.j*kx)
    return Hk
  
e1 = 2*np.pi/N
e2 = e1
kx = np.arange(0, 2*np.pi, e1)
ky = np.arange(0, 2*np.pi, e2)


S=0
for k1 in range(0, len(kx)):
    for k2 in range(0, len(ky)):

        w1, v1 = np.linalg.eig(HMat(phi, kx[k1], ky[k2]))
        idx1 = np.argsort(w1)
        v1_sorted = v1[:, idx1]

        w2, v2 = np.linalg.eig(HMat(phi, kx[k1]+e1, ky[k2]))
        idx2 = np.argsort(w2)
        v2_sorted = v2[:, idx2]

        w3, v3 = np.linalg.eig(HMat(phi, kx[k1], ky[k2]+e2))
        idx3 = np.argsort(w3)
        v3_sorted = v3[:, idx3]

        w4, v4 = np.linalg.eig(HMat(phi, kx[k1]+e1, ky[k2]+e2))
        idx4 = np.argsort(w4)
        v4_sorted = v4[:, idx4]

        U1 = np.linalg.det(v1_sorted)
        U1 = U1 / np.absolute(U1)
        U2 = np.linalg.det(v2_sorted)
        U2 = U2 / np.absolute(U2)
        U3 = np.linalg.det(v3_sorted)
        U3 = U3 / np.absolute(U3)
        U4 = np.linalg.det(v4_sorted) 
        U4 = U4 / np.absolute(U4)

        F = np.log(U1*U2*1/U3*1/U4)
        S = S+F

C = 1/(2*np.pi*1j)*S
C.real
