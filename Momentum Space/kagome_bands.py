import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt

t1 = -1
t2 = 0
#t2 = 0.3
#L1 = 0.28
L1 = -1
#L2 = 0.2
L2 = 0

def Hamiltonian(kx, ky):
    H = t1*np.array([
    [0, 1+np.exp(1j*kx), 1+np.exp(1j*ky)],
    [1+np.exp(-1j*kx), 0, 1+np.exp(-1j*(kx-ky))],
    [1+np.exp(-1j*ky), 1+np.exp(1j*(kx-ky)), 0]
    ])+1j*L1*np.array([
    [0, 1+np.exp(1j*kx), -(1+np.exp(1j*ky))],
    [-(1+np.exp(-1j*kx)), 0, 1+np.exp(-1j*(kx-ky))],
    [1+np.exp(-1j*ky), -(1+np.exp(1j*(kx-ky))), 0]
    ])+t2*np.array([
    [0, np.exp(1j*ky)+np.exp(1j*(kx-ky)), np.exp(1j*kx)+np.exp(-1j*(kx-ky))],
    [np.exp(-1j*ky)+np.exp(-1j*(kx-ky)), 0, np.exp(-1j*kx)+np.exp(-1j*ky)],
    [np.exp(-1j*kx)+np.exp(-1j*(kx-ky)), np.exp(1j*kx)+np.exp(-1j*ky), 0]
    ])+1j*L2*np.array([
    [0, -(np.exp(1j*ky)+np.exp(1j*(kx-ky))), np.exp(1j*kx)+np.exp(-1j*(kx-ky))],
    [np.exp(-1j*ky)+np.exp(-1j*(kx-ky)), 0, -(np.exp(-1j*kx)+np.exp(-1j*ky))],
    [-(np.exp(-1j*kx)+np.exp(-1j*(kx-ky))), np.exp(1j*kx)+np.exp(-1j*ky), 0]
    ])
    return H

size = len(Hamiltonian(1, 1))

def plot_dispersion(m,n):
    kx_range = np.linspace(np.pi, -np.pi, num=m)
    ky_range = np.linspace(np.pi, -np.pi, num=n)
    energies = np.zeros((m,n,size), dtype=complex)
    for i in range(m):
        for j in range(n):
            H = Hamiltonian(kx_range[i], ky_range[j])
            evals, evecs = LA.eigh(H)
            energies[j,i,:]=evals
    X, Y = np.meshgrid(kx_range, ky_range)
    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    ax.axis([X.min(), X.max(), Y.min(), Y.max()])
    for band in range(size):
        ax.plot_wireframe(X, Y, energies[:,:,band], rstride=1, cstride=1, color='0.5')
        #ax.scatter3D(X, Y, energies[:,:,band], antialiased=False, label="c")
    ax.set_xlabel(r'$k_{x}/\pi$', fontsize=15)
    ax.set_ylabel(r'$k_{y}/\pi$',fontsize=15)
    ax.set_zlabel(r'$E_k/|t_1|$', fontsize=15)
    ax.set_title('Momentum Uzayında Kagome Örgüsü \n'+r'$t_1=$'+str(t1)+','+r'$\lambda_1=$'+str(L1)+'\n'+r'$t_2=$'+str(t2)+','+r'$\lambda_2=$'+str(L2))
    #ax.view_init(azim=90, elev=0)
    ax.view_init(azim=45, elev=15)
    fig.savefig('kagome_momentum_3D_2.png', format='png', dpi=400)
    plt.show() 

plot_dispersion(25,25)
