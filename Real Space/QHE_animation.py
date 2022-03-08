#Libraries
from tight_binding_approximation import *
import matplotlib.animation as animation

#Hamiltonian Matrix of Hard-Wall Boundary Conditions
def HardHMat(alfa):
    H = np.zeros((L_x*L_y, L_x*L_y), dtype=complex)
    for m in range(L_x*L_y):
        for n in range(L_x*L_y):
            if m in HardBCLat[n]:
                if xy[m][0] > xy[n][0]:
                    H[m][n] = -np.exp(1j*2*np.pi*alfa*xy[m][1])
                elif xy[m][0] < xy[n][0]:
                    H[m][n] = -np.exp(-1j*2*np.pi*alfa*xy[m][1])
                else:
                    H[m][n]=-1
    return H

#Frame Per Second
fps = 5000 
#Frame Number of Animation
frn = 1500 

#Position of Inital Wave Function
start_lattice_size_1=5
start_lattice_size_2=6
#Normalization of State Vector
psi_0=np.zeros((L_x*L_y,1))/np.sqrt(2)
psi_0[start_lattice_size_1]=1
psi_0[start_lattice_size_2]=1

#Delta t
delta_t=0.01

#q Values Have To Be Eaxct Multiple of Lattice Size (Presence of Magnetic Field)
#alpha=0 (Absence of Magnetic Field)
alpha=1/(L_x/4)
#Series Expansion of Time Evolution Operator (Useful Way to Calculation)
time_evo=np.identity(L_x*L_y)-1.j*delta_t*HardHMat(alpha)

#Exact Time Evolution Operator (Memory Error!)
#time_evo_2=linalg.expm(-1.j*delta_t*HardHMat(alpha))

#Create Time Evolution Array
a = np.zeros((L_x*L_y,frn))
for i in range(1,frn):
    psi_0=np.matmul(time_evo, psi_0)
    a[:,[i]]=psi_0
    
#Create Animation
X = x_co
Y = y_co
xv, yv = np.meshgrid(X, Y)
def update_plot(frame_number, a, plot):
    plot[0].remove()
    plot[0] = ax.plot_surface(xv, yv, np.reshape((np.abs(a[:,[frame_number]]))**2, (L_x,L_y)), cmap="magma")
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.view_init(90,0)
plot = [ax.plot_surface(xv, yv, np.reshape((np.abs(a[:,[0]]))**2, (L_x,L_y)), color='0.75', rstride=1, cstride=1)]
ax.set_zlim(0,1.1)
ani = animation.FuncAnimation(fig, update_plot, frn, fargs=(a, plot), interval=1000/fps)

#Save Your Animation as .gif Format
#ani.save('animation.gif', writer='imagemagick', fps=200)
