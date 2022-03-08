#Libraries
from tight_binding_approximation import *
import sympy as sp
from matplotlib import cm
sp.init_printing(use_unicode=False, wrap_line=False, no_global=True)

#Magnetic Flux Per Plaquet
alpha=1/5

#Hamiltonian Matrix of Hard-Wall Boundary Conditions
def HardHMat(alpha):
    H = np.zeros((L_x*L_y, L_x*L_y), dtype=complex)
    for m in range(L_x*L_y):
        for n in range(L_x*L_y):
            if m in HardBCLat[n]:
                if xy[m][0] > xy[n][0]:
                    H[m][n] = -np.exp(1j*2*np.pi*alpha*xy[m][1])
                elif xy[m][0] < xy[n][0]:
                    H[m][n] = -np.exp(-1j*2*np.pi*alpha*xy[m][1])
                else:
                    H[m][n]=-1
    return H

#Symbolic Hamiltonian Matrix of Hard-Wall Boundary Conditions
def HardSHMat(alpha):
    H = sp.zeros(L_x*L_y, L_x*L_y, dtype=complex)
    for m in range(L_x*L_y):
        for n in range(L_x*L_y):
              if m in HardBCLat[n]:
                if xy[m,0] > xy[n,0]:
                    H[m,n] = -sp.exp(sp.I*2*sp.pi*alpha*xy[m][1])
                elif xy[m,0] < xy[n,0]:
                    H[m,n] = -sp.exp(-sp.I*2*sp.pi*alpha*xy[m][1])
                else:
                    H[m,n]=-1
    return H

#Hamiltonian Matrix of Periodic Boundary Conditions
def PerHMat(alpha):
    H = np.zeros((L_x*L_y, L_x*L_y), dtype=complex)
    for m in range(L_x*L_y):
        for n in range(L_x*L_y):
            if m in PerBCLat[n]:
                if np.absolute(xy[m][0]-xy[n][0])==L_x-1:
                    if xy[m][0] > xy[n][0]:
                        H[m][n] = -np.exp(-1j*2*np.pi*alpha*xy[m][1])
                    elif xy[m][0] < xy[n][0]:
                        H[m][n] = -np.exp(1j*2*np.pi*alpha*xy[m][1])
                else:
                    if xy[m][0] > xy[n][0]:
                        H[m][n] = -np.exp(1j*2*np.pi*alpha*xy[m][1])
                    elif xy[m][0] < xy[n][0]:
                        H[m][n] = -np.exp(-1j*2*np.pi*alpha*xy[m][1])
                    else:
                        H[m][n] = -np.exp(0)
    return H

#Symbolic Hamiltonian Matrix of Periodic Boundary Conditions
def PerSHMat(alpha):
    H = sp.zeros(L_x*L_y, L_x*L_y, dtype=complex)
    for m in range(L_x*L_y):
        for n in range(L_x*L_y):
              if m in PerBCLat[n]:
                if xy[m,0] > xy[n,0]:
                    H[m,n] = -sp.exp(sp.I*2*sp.pi*alpha*xy[m][1])
                elif xy[m,0] < xy[n,0]:
                    H[m,n] = -sp.exp(-sp.I*2*sp.pi*alpha*xy[m][1])
                else:
                    H[m,n]=-1
    return H

#Compare The Solutions of Hard-Wall and Periodic Boundary Conditions 
fig3, ax3 = plt.subplots()
y_real = np.zeros(L_x*L_y)
for i in range(L_x*L_y):
    y_real[i]=alpha
x_hard = np.linalg.eigvalsh(HardHMat(alpha))
x_per = np.linalg.eigvalsh(PerHMat(alpha))
plt.plot(x_hard, y_real, 'rx', markersize=7, label="Hard-Wall")
plt.plot(x_per, y_real, 'bo', markersize=7, label="Periodic")
plt.title('Solutions of Hard-Wall and Periodic B.C.')
plt.legend()
plt.show()

#Find Eigenvalues of Hamiltonian Matrix with Hard-Wall and Sorted It to Plot Sorted Eigenvalues
eigenValues, eigenVectors = np.linalg.eig(HardHMat(alpha))
idx = np.argsort(eigenValues)
sorted = eigenValues[idx]
fig5, ax5 = plt.subplots()
x = np.arange(0, L_x*L_y, 1)
y = sorted
plt.xlabel(r'$n$')
plt.ylabel(r'$E$')
plt.title('Sorted Eigenvalues of Hard-Wall Boundary Condition\n Size of lattice: '+str(L_x)+'x'+str(L_y))
plt.plot(x, y, 'ro', markersize=0.7)

# Plot Density of Nth. Eigenvector (You Can Find Edge States)
N=5
X = x_co
Y = y_co
xv, yv = np.meshgrid(X, Y)
Z = np.reshape((np.absolute(eigenVectors[:,idx[N]]))**2, (L_x,L_y))
fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
surf = ax.plot_surface(xv, yv, Z, cmap=cm.PRGn, linewidth=0, antialiased=False)
plt.xlabel('x')
plt.ylabel('y')
plt.title('Spesific Density of Choosen Eigenvector')
ax.set_zlabel(r'$|\Psi|^2$')
ax.view_init(50,30)
plt.colorbar(surf)

#Save Your Output
#plt.savefig("fig1.svg",format='svg', dpi=1200)
